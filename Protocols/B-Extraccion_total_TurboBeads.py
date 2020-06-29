import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv


# metadata
metadata = {
    'protocolName': 'Station B - RNA extraction',
    'author': 'Aitor Gastaminza & José Luis Villanueva & Alex Gasulla & Manuel Alba & Daniel Peñil',
    'source': 'Hospital Clínic Barcelona & HU Vall Hebrón & HU Marqués de Valdecilla',
    'apiLevel': '2.3',
    'description': 'Protocol for RNA extraction'
}

################################################
# CHANGE THESE VARIABLES ONLY
################################################
NUM_SAMPLES                         = 96    # Must be multiple of 8
LYSIS_VOLUME_PER_SAMPLE             = 300   # Original: 300
BEADS_VOLUME_PER_SAMPLE             = 420
WASH_VOLUME_PER_SAMPLE              = 300   # For each wash cycle
ELUTION_VOLUME_PER_SAMPLE           = 90
ELUTION_FINAL_VOLUME_PER_SAMPLE     = 50    # Volume transfered to final elution plate
LYSIS_NUM_MIXES                     = 5
BEADS_WELL_FIRST_TIME_NUM_MIXES     = 20
BEADS_WELL_NUM_MIXES                = 20
BEADS_NUM_MIXES                     = 20
WASH_NUM_MIXES                      = 20
ELUTION_NUM_MIXES                   = 20
VOLUME_SAMPLE                       = 200   # Sample volume received in station A
SET_TEMP_ON                         = True  # Do you want to start temperature module?
TEMPERATURE                         = 4     # Set temperature. It will be uesed if set_temp_on is set to True
################################################


run_id                      = 'B_Extraccion_total_TurboBeads'

recycle_tip                 = False # Do you want to recycle tips? It shoud only be set True for testing
mag_height                  = 7 # Height needed for NEST deepwell in magnetic deck
multi_well_rack_area        = 8 * 71 #Cross section of the 12 well reservoir

num_cols = math.ceil(NUM_SAMPLES / 8) # Columns we are working on

def run(ctx: protocol_api.ProtocolContext):

    #Change light to red
    ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))

    ctx.comment('Actual used columns: '+str(num_cols))
    STEP = 0
    STEPS = { #Dictionary with STEP activation, description, and times
            1:{'Execute': True, 'description': 'Transfer LYSIS'},   # Change the num of this step may affect remove supernatant step
            2:{'Execute': False, 'description': 'Wait rest', 'wait_time': 300},
            3:{'Execute': True, 'description': 'Transfer BEADS'},
            4:{'Execute': True, 'description': 'Wait rest', 'wait_time': 300},
            5:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 600}, 
            6:{'Execute': True, 'description': 'Remove supernatant'},
            7:{'Execute': True, 'description': 'Switch off magnet'},
            8:{'Execute': True, 'description': 'Add WASH'},
            9:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},
            10:{'Execute': True, 'description': 'Remove supernatant'},
            11:{'Execute': True, 'description': 'Switch off magnet'},
            12:{'Execute': True, 'description': 'Add WASH'},
            13:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},
            14:{'Execute': True, 'description': 'Remove supernatant'},
            15:{'Execute': True, 'description': 'Allow to dry', 'wait_time': 1200},
            16:{'Execute': True, 'description': 'Switch off magnet'},
            17:{'Execute': True, 'description': 'Add ELUTION'},
            18:{'Execute': False, 'description': 'Wait rest', 'wait_time': 300},
            19:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},
            20:{'Execute': True, 'description': 'Transfer to final elution plate'},
            }

    #Folder and file_path for log time
    import os
    folder_path = '/var/lib/jupyter/notebooks' + run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/Station_B_Extraccion_total_time_log.txt'

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
    Lysis = Reagent(name = 'Lysis',
                    flow_rate_aspirate = 1,
                    flow_rate_dispense = 1,
                    flow_rate_aspirate_mix = 1,
                    flow_rate_dispense_mix = 1,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = LYSIS_VOLUME_PER_SAMPLE, # reagent volume needed per sample
                    reagent_reservoir_volume =  (NUM_SAMPLES + 5) * LYSIS_VOLUME_PER_SAMPLE, 
                    num_wells = math.ceil((NUM_SAMPLES + 5) * LYSIS_VOLUME_PER_SAMPLE / 11500), #num_Wells max is 4, 13000 is the reservoir max volume (eventhough reservoir allows 15000)
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic
            
    Beads = Reagent(name = 'Beads',
                    flow_rate_aspirate = 3, #
                    flow_rate_dispense = 3, #
                    flow_rate_aspirate_mix = 25, #
                    flow_rate_dispense_mix = 50, #
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = BEADS_VOLUME_PER_SAMPLE, # reagent volume needed per sample
                    reagent_reservoir_volume =  math.ceil(NUM_SAMPLES * (BEADS_VOLUME_PER_SAMPLE + 100) * 1.1),  # 100 uL extra ispr per sample
                    num_wells = math.ceil(NUM_SAMPLES * (BEADS_VOLUME_PER_SAMPLE + 100) * 1.1 / 11500), #num_Wells max is 4, 13000 is the reservoir max volume (eventhough reservoir allows 15000)
                    h_cono = 1.95,
                    v_fondo = 695, #1.95 * multi_well_rack_area / 2, #Prismatic
                    tip_recycling = 'A1')

    Wash = Reagent(name = 'WASH',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 25,
                    flow_rate_dispense_mix = 100,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = WASH_VOLUME_PER_SAMPLE,
                    reagent_reservoir_volume = (2 * NUM_SAMPLES + 5) * WASH_VOLUME_PER_SAMPLE,
                    num_wells = math.ceil((2 * NUM_SAMPLES + 5) * WASH_VOLUME_PER_SAMPLE / 13000), 
                    h_cono = 1.95,
                    v_fondo = 695, #1.95 * multi_well_rack_area / 2, #Prismatic
                    tip_recycling = 'A1')

    Elution = Reagent(name = 'Elution',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 25,
                    flow_rate_dispense_mix = 40,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = False,
                    max_volume_allowed = 180,
                    reagent_volume = ELUTION_VOLUME_PER_SAMPLE,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * ELUTION_VOLUME_PER_SAMPLE,
                    num_wells = math.ceil((NUM_SAMPLES + 5) * ELUTION_VOLUME_PER_SAMPLE / 11500), #num_Wells max is 1
                    h_cono = 1.95,
                    v_fondo = 695) #1.95*multi_well_rack_area/2) #Prismatic

    Sample = Reagent(name = 'Sample',
                    flow_rate_aspirate = 0.5, # Original 0.5
                    flow_rate_dispense = 1, # Original 1
                    flow_rate_aspirate_mix = 1,
                    flow_rate_dispense_mix = 1,
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

    Lysis.vol_well      = Lysis.vol_well_original
    Beads.vol_well      = Beads.vol_well_original
    Wash.vol_well       = Wash.vol_well_original
    Elution.vol_well    = Elution.vol_well_original
    Sample.vol_well     = 350 # Arbitrary value

    #########
    def str_rounded(num):
        return str(int(num + 0.5))

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('VOLUMES FOR ' + str(NUM_SAMPLES) + ' SAMPLES')
    ctx.comment(' ')
    ctx.comment('Lysis: ' + str(Lysis.num_wells) + ' wells from well 2 in 12 well reservoir with volume ' + str_rounded(Lysis.vol_well_original) + ' uL each one')
    ctx.comment('Beads: ' + str(Beads.num_wells) + ' wells from well 6 in 12 well reservoir with volume ' + str_rounded(Beads.vol_well_original) + ' uL each one')
    ctx.comment('Elution: ' + str(Elution.num_wells) + ' wells from well 12 in 12 well reservoir with volume ' + str_rounded(Elution.vol_well_original) + ' uL each one')
    ctx.comment('Wash: in 195 mL reservoir 1 with volume ' + str(Wash.vol_well_original) + ' uL (+ dead volume)' )
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

    def move_vol_multi(pipet, reagent, source, dest, vol, x_offset_source, x_offset_dest, pickup_height, rinse, 
        avoid_droplet, wait_time, blow_out, touch_tip = False, touch_tip_v_offset = -10, drop_height = -5, 
        aspirate_with_x_scroll = False, dispense_bottom_air_gap_before = False):
        # Rinse before aspirating
        if rinse == True:
            custom_mix(pipet, reagent, location = source, vol = vol, rounds = 20, blow_out = False, mix_height = 3, offset = 0)

        # SOURCE
        if dispense_bottom_air_gap_before and reagent.air_gap_vol_bottom:
            pipet.dispense(reagent.air_gap_vol_bottom, source.top(z = -2), rate = reagent.flow_rate_dispense)

        if reagent.air_gap_vol_top != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_top) #air gap

        if aspirate_with_x_scroll:
            aspirate_with_x_scrolling(pip = pipet, volume = vol, src = source, pickup_height = pickup_height, rate = reagent.flow_rate_aspirate, start_x_offset_src = 0, stop_x_offset_src = x_offset_source)
        else:    
            s = source.bottom(pickup_height).move(Point(x = x_offset_source))
            pipet.aspirate(vol, s, rate = reagent.flow_rate_aspirate) # aspirate liquid

        if reagent.air_gap_vol_bottom != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_bottom) #air gap

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if avoid_droplet == True: # Touch the liquid surface to avoid droplets
            ctx.comment("Moving to: " + str(pickup_height))
            pipet.move_to(source.bottom(pickup_height))

        # GO TO DESTINATION
        d = dest.top(z = drop_height).move(Point(x = x_offset_dest))
        pipet.dispense(vol - reagent.disposal_volume + reagent.air_gap_vol_bottom, d, rate = reagent.flow_rate_dispense)

        if reagent.air_gap_vol_top != 0:
            pipet.dispense(reagent.air_gap_vol_top, dest.top(z = 0), rate = reagent.flow_rate_dispense)

        if blow_out == True:
            pipet.blow_out(dest.top(z = -5))

        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = touch_tip_v_offset, radius=0.7)
            
        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        #if reagent.air_gap_vol_bottom != 0:
            #pipet.move_to(dest.top(z = 0))
            #pipet.air_gap(reagent.air_gap_vol_bottom) #air gap
            #pipet.aspirate(air_gap_vol_bottom, dest.top(z = 0),rate = reagent.flow_rate_aspirate) #air gap

    def aspirate_with_x_scrolling(pip, volume, src, pickup_height = 0, rate = 1, start_x_offset_src = 0, stop_x_offset_src = 0):

        max_asp = volume/pip.min_volume
        inc_step = (start_x_offset_src - stop_x_offset_src) / max_asp

        for x in reversed(np.arange(stop_x_offset_src, start_x_offset_src, inc_step)):
            s = src.bottom(pickup_height).move(Point(x = x))
            pip.aspirate(volume = pip.min_volume, location = s, rate = rate)

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        #if not ctx.is_simulating():
        if tip_track['counts'][pip] >= tip_track['maxes'][pip]:
            for i in range(3):
                ctx._hw_manager.hardware.set_lights(rails=False)
                ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))
                time.sleep(0.3)
                ctx._hw_manager.hardware.set_lights(rails=True)
                ctx._hw_manager.hardware.set_lights(button=(0, 0 ,1))
                time.sleep(0.3)
            ctx._hw_manager.hardware.set_lights(button=(0, 1 ,0))
            ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
            resuming.')
            pip.reset_tipracks()
            tip_track['counts'][pip] = 0
            tip_track['num_refills'][pip] += 1
        pip.pick_up_tip()

    ##########
    def find_side(col):
        if col%2 == 0:
            side = -1 # left
        else:
            side = 1 # right
        return side

####################################
    # load labware and modules
    ######## 12 well rack
    reagent_res = ctx.load_labware('nest_12_reservoir_15ml', '7','reagent deepwell plate')

####################################
    ######## Single reservoirs
    reagent_res_1 = ctx.load_labware('nest_1_reservoir_195ml', '8', 'Single reagent reservoir 1')
    res_1 = reagent_res_1.wells()[0]

############################################
    ########## tempdeck
    tempdeck = ctx.load_module('Temperature Module Gen2', '1')
    if SET_TEMP_ON == True:
        tempdeck.set_temperature(TEMPERATURE)

##################################
    ####### Elution plate - final plate, goes to C
    #elution_plate = tempdeck.load_labware(
     #   'biorad_96_alum',
      #  'cooled elution plate')
    elution_plate = tempdeck.load_labware('kingfisher_96_aluminumblock_200ul', 
        'Kingfisher 96 Aluminum Block 200 uL')

############################################
    ######## Deepwell - comes from A
    magdeck = ctx.load_module('Magnetic Module Gen2', '4')
    #deepwell_plate = magdeck.load_labware('nest_96_wellplate_2ml_deep', 'NEST 96 Deepwell Plate 2mL') # Change to NEST deepwell plate.
    deepwell_plate = magdeck.load_labware('kingfisher_96_wellplate_2000ul', 'KingFisher 96 Well Plate 2mL') # Change to NEST deepwell plate.
    magdeck.disengage()

####################################
    ######## Waste reservoir
    waste_reservoir = ctx.load_labware('nest_1_reservoir_195ml', '11', 'waste reservoir') # Change to our waste reservoir
    waste = waste_reservoir.wells()[0] # referenced as reservoir

####################################
    ######### Load tip_racks
    tips300 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
        for slot in ['2', '3', '5', '6', '9']]

###############################################################################
    #Declare which reagents are in each reservoir as well as deepwell and elution plate
    Lysis.reagent_reservoir     = reagent_res.rows()[0][1:4]
    Beads.reagent_reservoir     = reagent_res.rows()[0][5:10]
    Elution.reagent_reservoir   = reagent_res.rows()[0][11:12]
    Wash.reagent_reservoir      = res_1
    work_destinations           = deepwell_plate.rows()[0][:Sample.num_wells]
    final_destinations          = elution_plate.rows()[0][:Sample.num_wells]

    # pipettes.
    m300 = ctx.load_instrument('p300_multi_gen2', 'right', tip_racks = tips300) # Load multi pipette

    #### used tip counter and set maximum tips available
    tip_track = {
        'counts': {m300: 0},
        'maxes': {m300: 96 * len(m300.tip_racks)}, #96 tips per tiprack * number or tipracks in the layout
        'num_refills' : {m300 : 0}
        }

###############################################################################

###############################################################################
    # STEP 1 TRANSFER LYSIS
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer lysis
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        lysis_trips = math.ceil(Lysis.reagent_volume / Lysis.max_volume_allowed)
        lysis_volume = Lysis.reagent_volume / lysis_trips
        lysis_transfer_vol = []
        for i in range(lysis_trips):
            lysis_transfer_vol.append(lysis_volume + Lysis.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False

        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j,transfer_vol in enumerate(lysis_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Lysis, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from reservoir column: ' + str(Lysis.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                move_vol_multi(m300, reagent = Lysis, source = Lysis.reagent_reservoir[Lysis.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True, touch_tip = False, drop_height = 1)
            
            if LYSIS_NUM_MIXES > 0:
                ctx.comment(' ')
                ctx.comment('Mixing sample ')
                custom_mix(m300, Lysis, location = work_destinations[i], vol =  Lysis.max_volume_allowed,
                        rounds = LYSIS_NUM_MIXES, blow_out = False, mix_height = 1, offset = 0)
            
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Lysis.air_gap_vol_bottom) #air gap
            
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
        # STEP 1 TRANSFER LYSIS
        ########

    ###############################################################################
    # STEP 2 WAIT REST
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Rest for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 2 WAIT REST
        ########

    ###############################################################################
    # STEP 3 TRANSFER BEADS
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer beads
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        beads_trips = math.ceil(Beads.reagent_volume / Beads.max_volume_allowed)
        beads_volume = Beads.reagent_volume / beads_trips
        beads_transfer_vol = []
        for i in range(beads_trips):
            beads_transfer_vol.append(beads_volume + Beads.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False # Original: True 
        first_mix_done = False

        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j,transfer_vol in enumerate(beads_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                transfer_vol_extra = transfer_vol if j > 0 else transfer_vol + 100  # Extra 100 isopropanol for calcs
                [pickup_height, change_col] = calc_height(Beads, multi_well_rack_area, transfer_vol_extra * 8)    
                if change_col == True or not first_mix_done: #If we switch column because there is not enough volume left in current reservoir column we mix new column
                    ctx.comment('Mixing new reservoir column: ' + str(Beads.col))
                    custom_mix(m300, Beads, Beads.reagent_reservoir[Beads.col],
                        vol = Beads.max_volume_allowed, rounds = BEADS_WELL_FIRST_TIME_NUM_MIXES, 
                        blow_out = False, mix_height = 1.5, offset = 0)
                    first_mix_done = True
                else:
                    ctx.comment('Mixing reservoir column: ' + str(Beads.col))
                    custom_mix(m300, Beads, Beads.reagent_reservoir[Beads.col],
                        vol = Beads.max_volume_allowed, rounds = BEADS_WELL_NUM_MIXES, 
                        blow_out = False, mix_height = 1.5, offset = 0)

                ctx.comment('Aspirate from reservoir column: ' + str(Beads.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                move_vol_multi(m300, reagent = Beads, source = Beads.reagent_reservoir[Beads.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True, touch_tip = True, drop_height = 1)
            
            if BEADS_NUM_MIXES > 0:
                ctx.comment(' ')
                ctx.comment('Mixing sample ')
                custom_mix(m300, Beads, location = work_destinations[i], vol =  Beads.max_volume_allowed,
                        rounds = BEADS_NUM_MIXES, blow_out = False, mix_height = 1, offset = 0, wait_time = 2)
            
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Beads.air_gap_vol_bottom) #air gap

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
        # STEP 3 TRANSFER BEADS
        ########

    ###############################################################################
    # STEP 4 WAIT REST
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Rest for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 4 WAIT REST
        ########

    ###############################################################################
    # STEP 5 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.comment(' ')
        magdeck.engage(height = mag_height)
        ctx.delay(seconds = STEPS[STEP]['wait_time'], msg = 'Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 5 INCUBATE WAIT WITH MAGNET ON
        ########

    ###############################################################################
    # STEP 6 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        actual_vol_well = Beads.reagent_volume + VOLUME_SAMPLE
        if STEPS[1]['Execute'] == True:             # Step 1 is lysis transfer
            actual_vol_well += Lysis.reagent_volume
        supernatant_trips = math.ceil((actual_vol_well) / Lysis.max_volume_allowed)
        supernatant_volume = Lysis.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Sample.disposal_volume)
        x_offset_rs = 2
        #Pickup_height is fixed here
        pickup_height = 0.5 # Original 0.5
        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            not_first_transfer = False

            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')

                move_vol_multi(m300, reagent = Sample, source = work_destinations[i],
                        dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = True,
                        dispense_bottom_air_gap_before = not_first_transfer)
                m300.air_gap(Sample.air_gap_vol_bottom)
                not_first_transfer = True

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
        # STEP 6 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 7 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 7 MAGNET OFF
        ########

    ###############################################################################
    # STEP 8 ADD WASH
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
        x_offset_rs = 2.5
        pickup_height = 0.5
        rinse = False # Not needed

        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in wash_transfer_vol:
                ctx.comment('Aspirate from reservoir 1')
                move_vol_multi(m300, reagent = Wash, source = Wash.reagent_reservoir,
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = False)
            
            if WASH_NUM_MIXES > 0:
                custom_mix(m300, Wash, location = work_destinations[i], vol = 180, two_thirds_mix_bottom = True,
                        rounds = WASH_NUM_MIXES, blow_out = False, mix_height = 3, offset = x_offset_dest)
            
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Wash.air_gap_vol_bottom) #air gap

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
        # STEP 8 ADD WASH
        ########

    ###############################################################################
    # STEP 9 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 9 INCUBATE WAIT WITH MAGNET ON
        ########

    ###############################################################################
    # STEP 10 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        supernatant_trips = math.ceil(Wash.reagent_volume / Wash.max_volume_allowed)
        supernatant_volume = Wash.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Sample.disposal_volume)
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            not_first_transfer = False

            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 0.5 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                move_vol_multi(m300, reagent = Sample, source = work_destinations[i],
                    dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                    pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False,
                    dispense_bottom_air_gap_before = not_first_transfer)
                m300.air_gap(Sample.air_gap_vol_bottom)
                not_first_transfer = True

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
        # STEP 10 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 11 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 11 MAGNET OFF
        ########

    ###############################################################################
    # STEP 12 ADD WASH
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
        x_offset_rs = 2.5
        pickup_height = 0.5
        rinse = False # Not needed

        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in wash_transfer_vol:
                ctx.comment('Aspirate from reservoir 1')
                move_vol_multi(m300, reagent = Wash, source = Wash.reagent_reservoir,
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = False)
            
            if WASH_NUM_MIXES > 0:
                custom_mix(m300, Wash, location = work_destinations[i], vol = 180, two_thirds_mix_bottom = True,
                    rounds = WASH_NUM_MIXES, blow_out = False, mix_height = 3, offset = x_offset_dest)
            
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Wash.air_gap_vol_bottom) #air gap

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
        # STEP 12 ADD WASH
        ########

    ###############################################################################
    # STEP 13 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 13 INCUBATE WAIT WITH MAGNET ON
        ########

    ###############################################################################
    # STEP 14 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        supernatant_trips = math.ceil(Wash.reagent_volume / Wash.max_volume_allowed)
        supernatant_volume = Wash.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Sample.disposal_volume)
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            not_first_transfer = False

            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 0.5 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                move_vol_multi(m300, reagent = Sample, source = work_destinations[i],
                    dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                    pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False,
                    dispense_bottom_air_gap_before = not_first_transfer)
                m300.air_gap(Sample.air_gap_vol_bottom)
                not_first_transfer = True

            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 14 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 15 ALLOW DRY
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')

        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Dry for ' + format(STEPS[STEP]['wait_time']) + ' seconds.') # 
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)
        ctx.comment('Used tips in total: ' + str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 15 ALLOW DRY
        ########


    ###############################################################################
    # STEP 16 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 16 MAGNET OFF
        ########
    
    ###############################################################################
    # STEP 17 ADD ELUTION
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
        elution_volume = Elution.reagent_volume / elution_trips
        elution_wash_vol = []
        for i in range(elution_trips):
            elution_wash_vol.append(elution_volume + Sample.disposal_volume)
        x_offset_rs = 2.5

        ########
        # Water or elution buffer
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs # Original 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in elution_wash_vol:
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Elution, multi_well_rack_area, transfer_vol*8)
                ctx.comment('Aspirate from Reservoir column: ' + str(Elution.col))
                ctx.comment('Pickup height is ' + str(pickup_height))

                move_vol_multi(m300, reagent = Elution, source = Elution.reagent_reservoir[Elution.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 0, blow_out = False, drop_height = -35)
            
            if ELUTION_NUM_MIXES > 0:
                ctx.comment(' ')
                ctx.comment('Mixing sample with Elution')
                custom_mix(m300, Elution, work_destinations[i], vol = Elution.reagent_volume, rounds = ELUTION_NUM_MIXES,
                    blow_out = False, mix_height = 1, offset = x_offset_dest, drop_height = -35)
            
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 17 ADD ELUTION
        ########

    ###############################################################################
    # STEP 18 WAIT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 18 WAIT
        ########

    ###############################################################################
    # STEP 19 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubate with magnet ON for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 19 INCUBATE WAIT WITH MAGNET ON
        ########

    ###############################################################################
    # STEP 20 TRANSFER TO ELUTION PLATE
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        elution_trips = math.ceil(ELUTION_FINAL_VOLUME_PER_SAMPLE / Elution.max_volume_allowed)
        elution_volume = ELUTION_FINAL_VOLUME_PER_SAMPLE / elution_trips
        elution_vol = []
        for i in range(elution_trips):
            elution_vol.append(elution_volume + Elution.disposal_volume)
        x_offset_rs = 2
        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in elution_vol:
                #Pickup_height is fixed here
                pickup_height = 1
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')

                move_vol_multi(m300, reagent = Sample, source = work_destinations[i],
                        dest = final_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = True, touch_tip = True)
            
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
                tip_track['counts'][m300] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 20 TRANSFER TO ELUTION PLATE
        ########

    '''if not ctx.is_simulating():
        with open(file_path,'w') as outfile:
            json.dump(STEPS, outfile)'''

    magdeck.disengage()
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

    # Light flash end of program
    import os
    #os.system('mpg123 /etc/audio/speaker-test.mp3')
    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        ctx._hw_manager.hardware.set_lights(button=(0, 0 ,1))
        time.sleep(0.3)
    ctx._hw_manager.hardware.set_lights(button=(0, 1 ,0))
    ctx.comment('Finished! \nMove deepwell plate (slot 5) to Station C for MMIX addition and PCR preparation.')
    used_tips = tip_track['num_refills'][m300] * 96 * len(m300.tip_racks) + tip_track['counts'][m300]
    ctx.comment('Used tips in total: '+str(used_tips))
    ctx.comment('Used racks in total: '+str(used_tips/96))
    ctx.comment('Available tips: '+str(tip_track['maxes'][m300]))
