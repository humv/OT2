import math
from opentrons.types import Point
from opentrons import protocol_api
import subprocess
import time
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'Station B - Kingfisher preparation',
    'author': 'Aitor Gastaminza & José Luis Villanueva & Alex Gasulla & Manuel Alba & Daniel Peñil',
    'source': 'Hospital Clínic Barcelona & HU Vall Hebrón & HU Marqués de Valdecilla',
    'apiLevel': '2.1',
    'description': 'Protocol for Kingfisher preparation'
}

################################################
# CHANGE THESE VARIABLES ONLY
################################################
NUM_SAMPLES                     = 96
BEADS_VOLUME_PER_SAMPLE         = 280
WASH_VOLUME_PER_SAMPLE          = 500
ETHANOL_VOLUME_PER_SAMPLE       = 500
ELUTION_VOLUME_PER_SAMPLE       = 50
BEADS_WELL_FIRST_TIME_NUM_MIXES = 10
BEADS_WELL_NUM_MIXES            = 3
BEADS_NUM_MIXES                 = 2

PHOTOSENSITIVE                  = False # True if it has photosensitive reagents
SOUND_NUM_PLAYS                 = 0
################################################

run_id                      = 'B-Preparacion_Kingfisher-Magmax_Viral_Pathogen'
path_sounds                 = '/var/lib/jupyter/notebooks/sonidos/'

recycle_tip     = False #
L_deepwell = 8 # Deepwell lenght (NEST deepwell)
#D_deepwell = 8.35 # Deepwell diameter (NUNC deepwell)
multi_well_rack_area = 8 * 71 #Cross section of the 12 well reservoir
deepwell_cross_section_area = L_deepwell ** 2 # deepwell square cross secion area

num_cols = math.ceil(NUM_SAMPLES / 8) # Columns we are working on
switch_off_lights = False

def run(ctx: protocol_api.ProtocolContext):

    #Change light to red
    ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))

    ctx.comment('Actual used columns: '+str(num_cols))
    STEP = 0
    STEPS = { #Dictionary with STEP activation, description, and times
            1:{'Execute': False, 'description': 'Transfer beads + PK + binding'},
            2:{'Execute': True, 'description': 'Transfer wash'},
            3:{'Execute': True, 'description': 'Transfer ethanol'},
            4:{'Execute': True, 'description': 'Transfer elution'}
            }

    #Folder and file_path for log time
    import os
    folder_path = '/var/lib/jupyter/notebooks/' + run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/Station_B_Preparacion_Kingfisher_time_log.txt'

    #Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, flow_rate_aspirate_mix, flow_rate_dispense_mix,
        air_gap_vol_bottom, air_gap_vol_top, disposal_volume, rinse, max_volume_allowed, reagent_volume, reagent_reservoir_volume, num_wells, h_cono, v_fondo, tip_recycling = 'none', dead_vol = 700):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
            self.flow_rate_dispense_mix = flow_rate_dispense_mix
            self.air_gap_vol_bottom = air_gap_vol_bottom
            self.air_gap_vol_top = air_gap_vol_top
            self.disposal_volume = disposal_volume
            self.rinse = bool(rinse)
            self.max_volume_allowed = max_volume_allowed
            self.reagent_volume = reagent_volume
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.num_wells = num_wells
            self.col = 0
            self.vol_well = 0
            self.h_cono = h_cono
            self.v_cono = v_fondo
            self.tip_recycling = tip_recycling
            self.dead_vol = dead_vol
            self.vol_well_original = (reagent_reservoir_volume / num_wells) + dead_vol if num_wells > 0 else 0

    #Reagents and their characteristics
    Wash = Reagent(name = 'Wash',
                    flow_rate_aspirate = 25, # Original = 0.5
                    flow_rate_dispense = 100, # Original = 1
                    flow_rate_aspirate_mix = 1, # Liquid density very high, needs slow aspiration
                    flow_rate_dispense_mix = 1, # Liquid density very high, needs slow dispensation
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = WASH_VOLUME_PER_SAMPLE, # reagent volume needed per sample
                    reagent_reservoir_volume =  (NUM_SAMPLES + 5) * WASH_VOLUME_PER_SAMPLE, #70000, #51648
                    num_wells = 1,
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    Ethanol = Reagent(name = 'Ethanol',
                    flow_rate_aspirate = 25,
                    flow_rate_dispense = 100,
                    flow_rate_aspirate_mix = 1,
                    flow_rate_dispense_mix = 1,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = ETHANOL_VOLUME_PER_SAMPLE,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * ETHANOL_VOLUME_PER_SAMPLE,
                    num_wells = 1, 
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    Beads_PK_Binding = Reagent(name = 'Magnetic beads + PK + Binding',
                    flow_rate_aspirate = 0.5,
                    flow_rate_dispense = 0.5,
                    flow_rate_aspirate_mix = 0.5,
                    flow_rate_dispense_mix = 0.5,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = BEADS_VOLUME_PER_SAMPLE,
                    reagent_reservoir_volume = NUM_SAMPLES * BEADS_VOLUME_PER_SAMPLE * 1.1,
                    num_wells = math.ceil(NUM_SAMPLES  * BEADS_VOLUME_PER_SAMPLE * 1.1 / 11500),
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    Elution = Reagent(name = 'Elution',
                    flow_rate_aspirate = 25, # Original 0.5
                    flow_rate_dispense = 100, # Original 1
                    flow_rate_aspirate_mix = 15,
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = False,
                    max_volume_allowed = 180,
                    reagent_volume = ELUTION_VOLUME_PER_SAMPLE,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * ELUTION_VOLUME_PER_SAMPLE,
                    num_wells = math.ceil((NUM_SAMPLES + 5) * ELUTION_VOLUME_PER_SAMPLE / 13000),
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    Sample = Reagent(name = 'Sample',
                    flow_rate_aspirate = 3, # Original 0.5
                    flow_rate_dispense = 3, # Original 1
                    flow_rate_aspirate_mix = 15,
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = False,
                    max_volume_allowed = 150,
                    reagent_volume = 50,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * 50, #14800,
                    num_wells = num_cols, #num_cols comes from available columns
                    h_cono = 4,
                    v_fondo = 4 * math.pi * 4**3 / 3) #Sphere

    Wash.vol_well               = Wash.vol_well_original
    Ethanol.vol_well            = Ethanol.vol_well_original
    Beads_PK_Binding.vol_well   = Beads_PK_Binding.vol_well_original
    Elution.vol_well            = Elution.vol_well_original
    Sample.vol_well             = 350 # Arbitrary value

    def str_rounded(num):
        return str(int(num + 0.5))

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('VOLUMES FOR ' + str(NUM_SAMPLES) + ' SAMPLES')
    ctx.comment(' ')
    if STEPS[1]['Execute']:
        ctx.comment('Beads + PK + Binding: ' + str(Beads_PK_Binding.num_wells) + ' wells from well 2 in multi reservoir with volume ' + str_rounded(Beads_PK_Binding.vol_well_original) + ' uL each one')
    ctx.comment('Elution: ' + str(Elution.num_wells) + ' wells from well 7 in multi reservoir with volume ' + str_rounded(Elution.vol_well_original) + ' uL each one')
    ctx.comment('Wash: in reservoir 1 with volume ' + str_rounded(Wash.vol_well_original) + ' uL')
    ctx.comment('Etanol: in reservoir 2 with volume ' + str_rounded(Ethanol.vol_well_original) + ' uL')
    ctx.comment('###############################################')
    ctx.comment(' ')

    ###################
    #Custom functions
    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height, offset, wait_time = 0, drop_height = -1, two_thirds_mix_bottom = False):
        '''
        Function for mix in the same location a certain number of rounds. Blow out optional. Offset
        can set to 0 or a higher/lower value which indicates the lateral movement
        '''
        if mix_height <= 0:
            mix_height = 1
        pipet.aspirate(1, location = location.bottom(z = mix_height), rate = reagent.flow_rate_aspirate_mix)
        for i in range(rounds):
            pipet.aspirate(vol, location = location.bottom(z = mix_height), rate = reagent.flow_rate_aspirate_mix)
            if two_thirds_mix_bottom and i < ((rounds / 3) * 2):
                pipet.dispense(vol, location = location.bottom(z = 5).move(Point(x = offset)), rate = reagent.flow_rate_dispense_mix)
            else:
                pipet.dispense(vol, location = location.top(z = drop_height).move(Point(x = offset)), rate = reagent.flow_rate_dispense_mix)
        pipet.dispense(1, location = location.bottom(z = mix_height), rate = reagent.flow_rate_dispense_mix)
        if blow_out == True:
            pipet.blow_out(location.top(z = -2)) # Blow out
        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

    def calc_height(reagent, cross_section_area, aspirate_volume, min_height = 0.4):
        nonlocal ctx
        ctx.comment('Remaining volume ' + str(reagent.vol_well) +
                    '< needed volume ' + str(aspirate_volume) + '?')
        if (reagent.vol_well - reagent.dead_vol) < aspirate_volume:
            ctx.comment('Next column should be picked')
            ctx.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            ctx.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            ctx.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
                    #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Remaining volume:' + str(reagent.vol_well))
            if height < min_height:
                height = min_height
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Calculated height is ' + str(height))
            if height < min_height:
                height = min_height
            ctx.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    def move_vol_multi(pipet, reagent, source, dest, vol, x_offset_source, x_offset_dest, pickup_height, rinse, avoid_droplet, wait_time, blow_out, touch_tip = False, drop_height = -5, dispense_bottom_air_gap_before = False):
        # Rinse before aspirating
        if rinse == True:
            #pipet.aspirate(air_gap_vol_top, location = source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap
            custom_mix(pipet, reagent, location = source, vol = vol, rounds = 20, blow_out = False, mix_height = 3, offset = 0)
            #pipet.dispense(air_gap_vol_top, location = source.top(z = -5), rate = reagent.flow_rate_dispense)

        # SOURCE
        if dispense_bottom_air_gap_before and reagent.air_gap_vol_bottom:
            pipet.dispense(reagent.air_gap_vol_bottom, source.top(z = -2), rate = reagent.flow_rate_dispense)


        if reagent.air_gap_vol_top != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_top) #air gap
            #pipet.aspirate(reagent.air_gap_vol_top, source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap

        s = source.bottom(pickup_height).move(Point(x = x_offset_source))
        pipet.aspirate(vol, s, rate = reagent.flow_rate_aspirate) # aspirate liquid

        if reagent.air_gap_vol_bottom != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.air_gap(reagent.air_gap_vol_bottom, height = 0) #air gap

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if avoid_droplet == True: # Touch the liquid surface to avoid droplets
            ctx.comment("Moving to: " + str(round(pickup_height, 2)) + ' mm')
            pipet.move_to(source.bottom(pickup_height))

        # GO TO DESTINATION
        d = dest.top(z = drop_height).move(Point(x = x_offset_dest))
        pipet.dispense(vol - reagent.disposal_volume + reagent.air_gap_vol_bottom, d, rate = reagent.flow_rate_dispense)

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if reagent.air_gap_vol_top != 0:
            pipet.dispense(reagent.air_gap_vol_top, dest.top(z = 0), rate = reagent.flow_rate_dispense)

        if blow_out == True:
            pipet.blow_out(dest.top(z = drop_height))

        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -10, radius=0.7)

        #if reagent.air_gap_vol_bottom != 0:
            #pipet.move_to(dest.top(z = 0))
            #pipet.air_gap(reagent.air_gap_vol_bottom) #air gap
            #pipet.aspirate(air_gap_vol_bottom, dest.top(z = 0),rate = reagent.flow_rate_aspirate) #air gap

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        #if not ctx.is_simulating():
        if tip_track['counts'][pip] >= tip_track['maxes'][pip]:
            ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
            resuming.')
            pip.reset_tipracks()
            tip_track['counts'][pip] = 0
        pip.pick_up_tip()

    ##########
    def start_run():
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Empezando protocolo')
        if PHOTOSENSITIVE == False:
            ctx._hw_manager.hardware.set_lights(button = True, rails =  True)
        else:
            ctx._hw_manager.hardware.set_lights(button = True, rails =  False)
        now = datetime.now()

        # dd/mm/YY H:M:S
        start_time = now.strftime("%Y/%m/%d %H:%M:%S")
        return start_time

    def run_quiet_process(command):
        subprocess.check_output('{} &> /dev/null'.format(command), shell=True)

    def play_sound(filename):
        print('Speaker')
        print('Next\t--> CTRL-C')
        try:
            run_quiet_process('mpg123 {}'.format(path_sounds + filename + '.mp3'))
        except KeyboardInterrupt:
            pass
            print()

    def finish_run(switch_off_lights = False):
        ctx.comment('###############################################')
        ctx.comment('Protocolo finalizado')
        ctx.comment(' ')
        #Set light color to blue
        ctx._hw_manager.hardware.set_lights(button = True, rails =  False)
        now = datetime.now()
        # dd/mm/YY H:M:S
        finish_time = now.strftime("%Y/%m/%d %H:%M:%S")
        if PHOTOSENSITIVE==False:
            for i in range(10):
                ctx._hw_manager.hardware.set_lights(button = False, rails =  False)
                time.sleep(0.3)
                ctx._hw_manager.hardware.set_lights(button = True, rails =  True)
                time.sleep(0.3)
        else:
            for i in range(10):
                ctx._hw_manager.hardware.set_lights(button = False, rails =  False)
                time.sleep(0.3)
                ctx._hw_manager.hardware.set_lights(button = True, rails =  False)
                time.sleep(0.3)
        if switch_off_lights:
            ctx._hw_manager.hardware.set_lights(button = True, rails =  False)

        # TODO: Añadir refills a los tip_racks
        # used_tips = tip_track['num_refills'][m300] * 96 * len(m300.tip_racks) + tip_track['counts'][m300]
        ctx.comment('Puntas de 200 uL utilizadas: ' + str(tip_track['counts'][m300]) + ' (' + str(round(tip_track['counts'][m300] / 96, 2)) + ' caja(s))')
        ctx.comment('###############################################')

        if not ctx.is_simulating():
            for i in range(SOUND_NUM_PLAYS):
                if i > 0:
                    time.sleep(60)
                play_sound('finalizado')

        return finish_time


    ##########
    def find_side(col):
        if col%2 == 0:
            side = -1 # left
        else:
            side = 1 # right
        return side

# load labware and modules
############################################
    ######## Sample plate - comes from A
    deepwell_plate_samples = ctx.load_labware('nest_96_wellplate_2ml_deep', '1' , 'NEST 96 Deepwell Plate 2mL 1')
    deepwell_plate_wash = ctx.load_labware('nest_96_wellplate_2ml_deep', '5' , 'NEST 96 Deepwell Plate 2mL 2') 
    deepwell_plate_ethanol = ctx.load_labware('nest_96_wellplate_2ml_deep', '6' , 'NEST 96 Deepwell Plate 2mL 3') 
    deepwell_plate_elution = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '7' , 'Bio-Rad 96 Well Plate 200 µL PCR')

####################################
    ######## 12 well rack
    reagent_multi_res = ctx.load_labware('nest_12_reservoir_15ml', '4','Reagent deepwell plate')

####################################
    ######## Single reservoirs
    reagent_res_1 = ctx.load_labware('nest_1_reservoir_195ml', '2', 'Single reagent reservoir 1')
    res_1 = reagent_res_1.wells()[0]

    reagent_res_2 = ctx.load_labware('nest_1_reservoir_195ml', '3', 'Single reagent reservoir 2')
    res_2 = reagent_res_2.wells()[0]

####################################
    ######### Load tip_racks
    tips300 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
        for slot in ['8', '9']]

###############################################################################
    #Declare which reagents are in each reservoir as well as deepwell and sample plate
    Wash.reagent_reservoir      = res_1
    Ethanol.reagent_reservoir   = res_2
    Beads_PK_Binding.reagent_reservoir  = reagent_multi_res.rows()[0][1:5]
    Elution.reagent_reservoir   = reagent_multi_res.rows()[0][6:7]
    work_destinations           = deepwell_plate_samples.rows()[0][:Sample.num_wells]
    wash_destinations           = deepwell_plate_wash.rows()[0][:Sample.num_wells]
    ethanol_destinations        = deepwell_plate_ethanol.rows()[0][:Sample.num_wells]
    elution_destinations        = deepwell_plate_elution.rows()[0][:Sample.num_wells]

    # pipettes.
    m300 = ctx.load_instrument('p300_multi_gen2', 'right', tip_racks=tips300) # Load multi pipette

    #### used tip counter and set maximum tips available
    tip_track = {
        'counts': {m300: 0},
        'maxes': {m300: 96 * len(m300.tip_racks)} #96 tips per tiprack * number or tipracks in the layout
        }

###############################################################################

###############################################################################
    start_run()
    
    ###############################################################################
    # STEP 1 TRANSFER BEADS + PK + Binding
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        beads_trips = math.ceil(Beads_PK_Binding.reagent_volume / Beads_PK_Binding.max_volume_allowed)
        beads_volume = Beads_PK_Binding.reagent_volume / beads_trips #136.66
        beads_transfer_vol = []
        for i in range(beads_trips):
            beads_transfer_vol.append(beads_volume + Beads_PK_Binding.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False # Original: True
        first_mix_done = False

        for i in range(num_cols):
            not_first_transfer = False

            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j,transfer_vol in enumerate(beads_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Beads_PK_Binding, multi_well_rack_area, transfer_vol * 8)
                if change_col == True or not first_mix_done: #If we switch column because there is not enough volume left in current reservoir column we mix new column
                    ctx.comment('Mixing new reservoir column: ' + str(Beads_PK_Binding.col))
                    custom_mix(m300, Beads_PK_Binding, Beads_PK_Binding.reagent_reservoir[Beads_PK_Binding.col],
                            vol = Beads_PK_Binding.max_volume_allowed, rounds = BEADS_WELL_FIRST_TIME_NUM_MIXES, blow_out = False, mix_height = 0.5, offset = 0)
                    first_mix_done = True
                else:
                    ctx.comment('Mixing reservoir column: ' + str(Beads_PK_Binding.col))
                    custom_mix(m300, Beads_PK_Binding, Beads_PK_Binding.reagent_reservoir[Beads_PK_Binding.col],
                            vol = Beads_PK_Binding.max_volume_allowed, rounds = BEADS_WELL_NUM_MIXES, blow_out = False, mix_height = 0.5, offset = 0)
                ctx.comment('Aspirate from reservoir column: ' + str(Beads_PK_Binding.col))
                ctx.comment('Pickup height is ' + str(round(pickup_height, 2)) + ' mm')
                #if j!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = Beads_PK_Binding, source = Beads_PK_Binding.reagent_reservoir[Beads_PK_Binding.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 2, blow_out = False, 
                        touch_tip = False, drop_height = 5, dispense_bottom_air_gap_before = not_first_transfer)
                
                m300.air_gap(Beads_PK_Binding.air_gap_vol_bottom, height = 5)
                not_first_transfer = True
            
            ctx.comment(' ')
            ctx.comment('Mixing sample ')
            custom_mix(m300, Beads_PK_Binding, location = work_destinations[i], vol =  Beads_PK_Binding.max_volume_allowed,
                    rounds = BEADS_NUM_MIXES, blow_out = False, mix_height = 0, offset = 0, wait_time = 2, two_thirds_mix_bottom = True)
            m300.air_gap(Beads_PK_Binding.air_gap_vol_bottom, height = 0) #air gap

            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 1 TRANSFER BEADS + PK + Binding
        ########

    ###############################################################################
    # STEP 2 TRANSFER WASH
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        wash_trips = math.ceil(Wash.reagent_volume / Wash.max_volume_allowed)
        wash_volume = Wash.reagent_volume / wash_trips #136.66
        wash_transfer_vol = []
        for i in range(wash_trips):
            wash_transfer_vol.append(wash_volume + Wash.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False
        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)

            for j,transfer_vol in enumerate(wash_transfer_vol):
                ctx.comment('Aspirate from reservoir 1')
                #if j!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = Wash, source = Wash.reagent_reservoir,
                        dest = wash_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = 2, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True, touch_tip = False)
            ctx.comment(' ')
        if recycle_tip == True:
            m300.return_tip()
        else:
            m300.drop_tip(home_after = False)
        tip_track['counts'][m300] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 2 TRANSFER WASH
        ########

    ###############################################################################
    # STEP 3 TRANSFER ETHANOL
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ethanol_trips = math.ceil(Ethanol.reagent_volume / Ethanol.max_volume_allowed)
        ethanol_volume = Ethanol.reagent_volume / ethanol_trips #136.66
        ethanol_transfer_vol = []
        for i in range(ethanol_trips):
            ethanol_transfer_vol.append(ethanol_volume + Ethanol.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False
        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)

            for j,transfer_vol in enumerate(ethanol_transfer_vol):
                ctx.comment('Aspirate from reservoir 2')
                #if j!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = Ethanol, source = Ethanol.reagent_reservoir,
                        dest = ethanol_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = 2, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True, touch_tip = False)

        if recycle_tip == True:
            m300.return_tip()
        else:
            m300.drop_tip(home_after = False)
        tip_track['counts'][m300] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 3 TRANSFER ETHANOL
        ########

    ###############################################################################
    # STEP 4 TRANSFER ELUTION
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        elution_trips = math.ceil(Elution.reagent_volume / Elution.max_volume_allowed)
        elution_volume = Elution.reagent_volume / elution_trips #136.66
        elution_transfer_vol = []
        for i in range(elution_trips):
            elution_transfer_vol.append(elution_volume + Elution.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False

        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j,transfer_vol in enumerate(elution_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Elution, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from reservoir column: ' + str(Elution.col))
                ctx.comment('Pickup height is ' + str(round(pickup_height, 2)) + ' mm')
                #if j!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = Elution, source = Elution.reagent_reservoir[Elution.col],
                        dest = elution_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = False, touch_tip = False)
            ctx.comment(' ')

        if recycle_tip == True:
            m300.return_tip()
        else:
            m300.drop_tip(home_after = False)
        tip_track['counts'][m300] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 4 TRANSFER ELUTION
        ########


    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('Homing robot')
    ctx.comment('###############################################')
    ctx.comment(' ')
    ctx.home()
###############################################################################
    # Export the time log to a tsv file
    if not ctx.is_simulating():
        with open(file_path, 'w') as f:
            f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
            for key in STEPS.keys():
                row = str(key)
                for key2 in STEPS[key].keys():
                    row += '\t' + format(STEPS[key][key2])
                f.write(row + '\n')
        f.close()

    finish_run(switch_off_lights)

