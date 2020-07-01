import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv


# metadata
metadata = {
    'protocolName': 'Station B - RNA extraction',
    'author': 'Aitor Gastaminza & José Luis Villanueva & Alex Gasulla & Manuel Alba & Daniel Peñil & Francisco Martinez (HU Virgen del Rocio) & Rafael Martinez (HU Virgen del Rocio)',
    'source': 'Hospital Clínic Barcelona & HU Vall Hebrón & HU Marqués de Valdecilla & HU Virgen del Rocio',
    'apiLevel': '2.4',
    'description': 'Protocol for RNA extraction'
}

################################################
# CHANGE THESE VARIABLES ONLY
################################################
NUM_SAMPLES                 = 8
FAGO_VOLUME_PER_SAMPLE      = 20
LYSIS_VOLUME_PER_SAMPLE     = 280
WASH_VOLUME_PER_SAMPLE      = 500
ETHANOL_VOLUME_PER_SAMPLE   = 500
ELUTION_VOLUME_PER_SAMPLE   = 50
VOLUME_SAMPLE               = 300   # Sample volume received in station A
SET_TEMP_ON                 = False  # Do you want to start temperature module?
TEMPERATURE                 = 4     # Set temperature. It will be uesed if set_temp_on is set to True
FIRST_TIPS_COLUMN           = 0
WAIT_MAGNET                 = 600
################################################

RECYCLE_TIP                 = False # Do you want to recycle tips? It shoud only be set True for testing

run_id                      = 'B_Extraccion_total'

#mag_height = 11 # Height needed for NUNC deepwell in magnetic deck
mag_height                  = 7 # Height needed for NEST deepwell in magnetic deck

L_deepwell                  = 7.5 # Deepwell lenght (NEST deepwell)
#D_deepwell = 8.35 # Deepwell diameter (NUNC deepwell)
multi_well_rack_area        = 8 * 71 #Cross section of the 12 well reservoir
deepwell_cross_section_area = L_deepwell ** 2 # deepwell square cross secion area

num_cols = math.ceil(NUM_SAMPLES / 8) # Columns we are working on

def run(ctx: protocol_api.ProtocolContext):

    #Change light to red
    ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))

    ctx.comment('Actual used columns: ' + str(num_cols))
    STEP = 0
    STEPS = { #Dictionary with STEP activation, description, and times
         1:{'Execute': True, 'description': 'Transfer FAGO'},#
         2:{'Execute': True, 'description': 'Transfer LYSIS'},#
         3:{'Execute': False, 'description': 'Wait rest', 'wait_time': 300},# 300
         4:{'Execute': False, 'description': 'Incubate wait with magnet ON', 'wait_time': WAIT_MAGNET}, #600
         5:{'Execute': True, 'description': 'Remove supernatant'},#
         6:{'Execute': False, 'description': 'Switch off magnet'},#
         7:{'Execute': True, 'description': 'Add WASH'},#
         8:{'Execute': False, 'description': 'Incubate wait with magnet ON', 'wait_time': WAIT_MAGNET},#300
         9:{'Execute': True, 'description': 'Remove supernatant'},#
        10:{'Execute': False, 'description': 'Switch off magnet'},#
        11:{'Execute': True, 'description': 'Add ETHANOL'},#
        12:{'Execute': False, 'description': 'Incubate wait with magnet ON', 'wait_time': WAIT_MAGNET},#300
        13:{'Execute': True, 'description': 'Remove supernatant'},#
        14:{'Execute': False, 'description': 'Allow to dry', 'wait_time': WAIT_MAGNET},#600
        15:{'Execute': False, 'description': 'Switch off magnet'},#
        16:{'Execute': True, 'description': 'Add ELUTION'},#
        17:{'Execute': False, 'description': 'Wait rest', 'wait_time': 600},#60
        18:{'Execute': False, 'description': 'Incubate wait with magnet ON', 'wait_time': WAIT_MAGNET},#300
        19:{'Execute': True, 'description': 'Transfer to final elution plate'},
    }

    #Folder and file_path for log time
    import os
    folder_path = '/var/lib/jupyter/notebooks' + run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/Station_B_Extraccion_total_time_log.txt'

    #Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, flow_rate_aspirate_mix, flow_rate_dispense_mix,
        air_gap_vol_bottom, air_gap_vol_top, disposal_volume, rinse, max_volume_allowed, reagent_volume, reagent_reservoir_volume, num_wells, h_cono, v_fondo, tip_recycling = 'none'):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
            self.flow_rate_dispense_mix = flow_rate_dispense_mix
            self.air_gap_vol_bottom = air_gap_vol_bottom
            self.air_gap_vol_top = air_gap_vol_top
            self.disposal_volume = disposal_volume
            self.rinse = bool(rinse)
            self.max_volume_allowed = max_volume_allowed
            self.reagent_volume = reagent_volume
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.num_wells = num_wells
            self.col = 0
            self.vol_well = 0
            self.h_cono = h_cono
            self.v_cono = v_fondo
            self.tip_recycling = tip_recycling
            self.vol_well_original = reagent_reservoir_volume / num_wells

    #Reagents and their characteristics
    Fago = Reagent(name='Internal Control',
        flow_rate_aspirate=1,
        flow_rate_dispense=1,
        flow_rate_aspirate_mix = 1,
        flow_rate_dispense_mix = 1,
        air_gap_vol_bottom = 5,
        air_gap_vol_top = 5,
        disposal_volume = 1,
        rinse=True,
        max_volume_allowed = 180,
        reagent_volume = FAGO_VOLUME_PER_SAMPLE, # reagent volume needed per sample
        reagent_reservoir_volume =  (NUM_SAMPLES + 5) * FAGO_VOLUME_PER_SAMPLE, 
        num_wells = 1,
        h_cono = 1.95,
        v_fondo = 750,
        tip_recycling = 'A1')

    Lysis = Reagent(name = 'Lysis',
        flow_rate_aspirate = 0.5, # Original = 0.5
        flow_rate_dispense = 0.5, # Original = 1
        flow_rate_aspirate_mix = 0.5, # Liquid density very high, needs slow aspiration
        flow_rate_dispense_mix = 0.5, # Liquid density very high, needs slow dispensation
        air_gap_vol_bottom = 5,
        air_gap_vol_top = 5,
        disposal_volume = 1,
        rinse = True,
        max_volume_allowed = 180,
        reagent_volume = LYSIS_VOLUME_PER_SAMPLE, # reagent volume needed per sample
        reagent_reservoir_volume =  (NUM_SAMPLES + 3) * LYSIS_VOLUME_PER_SAMPLE, 
        num_wells = math.ceil((NUM_SAMPLES + 3) * LYSIS_VOLUME_PER_SAMPLE / 10300), #num_Wells max is 4, 13000 is the reservoir max volume (eventhough reservoir allows 15000)
        h_cono = 1.95,
        v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
        tip_recycling = 'A1')

    Wash = Reagent(name = 'WASH',
        flow_rate_aspirate = 3,
        flow_rate_dispense = 3, #3
        flow_rate_aspirate_mix = 15,
        flow_rate_dispense_mix = 25,
        air_gap_vol_bottom = 5,#5
        air_gap_vol_top = 5, #5
        disposal_volume = 1,
        rinse = True,
        max_volume_allowed = 180,
        reagent_volume = WASH_VOLUME_PER_SAMPLE,
        reagent_reservoir_volume = (NUM_SAMPLES + 5) * WASH_VOLUME_PER_SAMPLE,
        num_wells = math.ceil((NUM_SAMPLES + 5) * WASH_VOLUME_PER_SAMPLE / 10300), 
        h_cono = 1.95,
        v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
        tip_recycling = 'A1')

    Ethanol = Reagent(name = 'Ethanol',
        flow_rate_aspirate = 1,
        flow_rate_dispense = 1,
        flow_rate_aspirate_mix = 5,
        flow_rate_dispense_mix = 5,
        air_gap_vol_bottom = 5,
        air_gap_vol_top = 5,
        disposal_volume = 1,
        rinse = True,
        max_volume_allowed = 180,
        reagent_volume = ETHANOL_VOLUME_PER_SAMPLE,
        reagent_reservoir_volume = (NUM_SAMPLES + 5) * ETHANOL_VOLUME_PER_SAMPLE,
        num_wells = math.ceil((NUM_SAMPLES + 5) * ETHANOL_VOLUME_PER_SAMPLE / 10300), 
        h_cono = 1.95,
        v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
        tip_recycling = 'A1')

    Elution = Reagent(name = 'Elution',
        flow_rate_aspirate = 1,
        flow_rate_dispense = 1,
        flow_rate_aspirate_mix = 0.5,
        flow_rate_dispense_mix = 0.5,
        air_gap_vol_bottom = 5,
        air_gap_vol_top = 5,
        disposal_volume = 1,
        rinse = False,
        max_volume_allowed = 180,
        reagent_volume = 50,
        reagent_reservoir_volume = (NUM_SAMPLES + 5) * ELUTION_VOLUME_PER_SAMPLE,
        num_wells = math.ceil((NUM_SAMPLES + 5) * ELUTION_VOLUME_PER_SAMPLE / 10300),
        h_cono = 1.95,
        v_fondo = 750) #1.95*multi_well_rack_area/2) #Prismatic

    Sample = Reagent(name = 'Sample',
        flow_rate_aspirate = 3, # Original 0.5
        flow_rate_dispense = 3, # Original 1
        flow_rate_aspirate_mix = 15,
        flow_rate_dispense_mix = 25,
        air_gap_vol_bottom = 5,
        air_gap_vol_top = 5,
        disposal_volume = 1,
        rinse = False,
        max_volume_allowed = 150,
        reagent_volume = 50,
        reagent_reservoir_volume = (NUM_SAMPLES + 5) * 50, #14800,
        num_wells = num_cols, #num_cols comes from available columns
        h_cono = 4,
        v_fondo = 4 * math.pi * 4**3 / 3) #Sphere

    Fago.vol_well       = Fago.vol_well_original
    Lysis.vol_well      = Lysis.vol_well_original
    Wash.vol_well       = Wash.vol_well_original
    Ethanol.vol_well    = Ethanol.vol_well_original
    Elution.vol_well    = Elution.vol_well_original
    Sample.vol_well     = 300 # Arbitrary value

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('VOLUMES FOR ' + str(NUM_SAMPLES) + ' SAMPLES')
    ctx.comment(' ')
    ctx.comment('Fago: ' + str(Fago.num_wells) + ' well 2 in 2 well reservoir with volume ' + str(Fago.vol_well_original) + ' uL each one')
    ctx.comment('Lysis: ' + str(Lysis.num_wells) + ' wells from 4 to 8 in 12 well reservoir with volume ' + str(Lysis.vol_well_original) + ' uL each one')
    ctx.comment('Elution: ' + str(Elution.num_wells) + ' wells from 10 to 11 in 12 well reservoir with volume ' + str(Elution.vol_well_original) + ' uL each one')
    ctx.comment('Wash: in 195 mL reservoir 1 with volume ' + str(Wash.vol_well_original) + ' uL')
    ctx.comment('Ethanol: in 195 mL reservoir 2 with volume ' + str(Ethanol.vol_well_original) + ' uL')
    ctx.comment('###############################################')
    ctx.comment(' ')

    ###################
    #Custom functions
    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height, offset, wait_time = 0, mix_height_bot = -5):
        '''
        Function for mix in the same location a certain number of rounds. Blow out optional. Offset
        can set to 0 or a higher/lower value which indicates the lateral movement
        '''
        if mix_height == 0:
            mix_height = 1
        pipet.aspirate(1, location = location.bottom(z = mix_height), rate = reagent.flow_rate_aspirate_mix)
        for _ in range(rounds):
            pipet.aspirate(vol, location = location.bottom(z = mix_height), rate = reagent.flow_rate_aspirate_mix)
            pipet.dispense(vol, location = location.top(z = mix_height_bot).move(Point(x = offset)), rate = reagent.flow_rate_dispense_mix)
        pipet.dispense(1, location = location.bottom(z = mix_height), rate = reagent.flow_rate_dispense_mix)
        if blow_out == True:
            pipet.blow_out(location.top(z = -2)) # Blow out
        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

    def calc_height(reagent, cross_section_area, aspirate_volume):
        nonlocal ctx
        ctx.comment('Remaining volume ' + str(reagent.vol_well) +
                    '< needed volume ' + str(aspirate_volume) + '?')
        if reagent.vol_well < aspirate_volume:
            ctx.comment('Next column should be picked')
            ctx.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            ctx.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            ctx.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
                    #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Remaining volume:' + str(reagent.vol_well))
            if height < 5:
                height = 1
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Calculated height is ' + str(height))
            if height < 5:
                height = 1
            ctx.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    def move_vol_multi(pipet, reagent, source, dest, vol, x_offset_source, x_offset_dest, pickup_height, rinse,
        avoid_droplet, wait_time, blow_out, touch_tip = False, drop_height = -5, blow_height = -5, blow_wash=False):
        # Rinse before aspirating
        if rinse == True:
            #pipet.aspirate(air_gap_vol_top, location = source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap
            custom_mix(pipet, reagent, location = source, vol = vol, rounds = 20, blow_out = False, mix_height = 3, offset = 0)
            #pipet.dispense(air_gap_vol_top, location = source.top(z = -5), rate = reagent.flow_rate_dispense)

        # SOURCE
        if reagent.air_gap_vol_top != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_top) #air gap
            #pipet.aspirate(reagent.air_gap_vol_top, source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap

        s = source.bottom(pickup_height).move(Point(x = x_offset_source))
        pipet.aspirate(vol, s) # aspirate liquid

        if reagent.air_gap_vol_bottom != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_bottom) #air gap
            #pipet.aspirate(air_gap_vol_bottom, source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if avoid_droplet == True: # Touch the liquid surface to avoid droplets
            ctx.comment("Moving to: " + str(pickup_height))
            pipet.move_to(source.bottom(pickup_height))

        # GO TO DESTINATION
        d = dest.top(z = drop_height).move(Point(x = x_offset_dest))
        pipet.dispense(vol - reagent.disposal_volume + reagent.air_gap_vol_bottom, d, rate = reagent.flow_rate_dispense)

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if reagent.air_gap_vol_top != 0:
            pipet.dispense(reagent.air_gap_vol_top, dest.top(z = 0), rate = reagent.flow_rate_dispense)

        if blow_out == True:
            pipet.blow_out(dest.top(z = blow_height))

        if blow_wash == True:
            pipet.air_gap(150)

        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -10, radius=0.7)

        #if reagent.air_gap_vol_bottom != 0:
            #pipet.move_to(dest.top(z = 0))
            #pipet.air_gap(reagent.air_gap_vol_bottom) #air gap
            #pipet.aspirate(air_gap_vol_bottom, dest.top(z = 0),rate = reagent.flow_rate_aspirate) #air gap

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip, tips):
        nonlocal tip_track
        if tip_track['actual'][tips] >= tip_track['maxes'][tips]:
            ctx.pause('Replace ' + str(tips) + ' before resuming.')
            pip.reset_tipracks()
            tip_track['actual'][tips] = 0
            tip_track['maxes'][tips] = 96
            FIRST_TIPS_COLUMN = 0
        pip.pick_up_tip(tips.wells()[tip_track['actual'][tips]])
        tip_track['actual'][tips] += 8 

    def drop_tip(pip, tips, drop):
        if (RECYCLE_TIP):
            pip.return_tip(home_after = False)
            tip_track['counts'][tips] += 8
        elif(drop):
            pip.drop_tip(home_after = False)
            tip_track['counts'][tips] += 8
        else:
            pip.return_tip(home_after = False)

    ##########
    def find_side(col):
        if col%2 == 0:
            side = -0.9 # left
        else:
            side = 0.9 # right
        return side

####################################
    # load labware and modules
    ######## 12 well rack
    reagent_res_1 = ctx.load_labware('nest_12_reservoir_15ml', '2','Agentes (FAGO / Lysis / Elution)')
    reagent_res_2 = ctx.load_labware('nest_12_reservoir_15ml', '1','Agentes (Wash / Ethanol)')

############################################
    ########## tempdeck
    tempdeck = ctx.load_module('Temperature Module Gen2', '3')

##################################
    ####### Elution plate - final plate, goes to C
    elution_plate = tempdeck.load_labware('opentrons_96_aluminumblock_nest_wellplate_100ul', 'PCR para C')

############################################
    ######## Elution plate - comes from A
    magdeck = ctx.load_module('Magnetic Module Gen2', '6')
    deepwell_plate = magdeck.load_labware('thermo_96_wellplate_2200ul', 'DeepWeel de A') # Change to NEST deepwell plate.
    # deepwell_plate = magdeck.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 'thermo_96_wellplate_2200ul') # Change to NEST deepwell plate.
    magdeck.disengage()

####################################
    ######## Waste reservoir
    waste_reservoir = ctx.load_labware('nest_1_reservoir_195ml', '9', 'Liquidos sobrantes') # Change to our waste reservoir
    waste = waste_reservoir.wells()[0] # referenced as reservoir

####################################
    ######### Load tip_racks
    tips300Beads     = ctx.load_labware('opentrons_96_tiprack_300ul', '4', 'Puntas 200µl para Magnet')
    tips300Fago      = ctx.load_labware('opentrons_96_tiprack_300ul', '5', 'Puntas 200µl para Fago')
    tips300Wash      = ctx.load_labware('opentrons_96_tiprack_300ul', '7', 'Puntas 200µl para Wash')
    tips300Ethanol   = ctx.load_labware('opentrons_96_tiprack_300ul', '8', 'Puntas 200µl para Ethanol')
    tips300Elution   = ctx.load_labware('opentrons_96_tiprack_300ul', '10', 'Puntas 200µl para Elution')
    tips300Samples   = ctx.load_labware('opentrons_96_tiprack_300ul', '11', 'Puntas 200µl para Samples')

###############################################################################
    #Declare which reagents are in each reservoir as well as deepwell and elution plate
    Fago.reagent_reservoir       = reagent_res_1.rows()[0][1]
    Lysis.reagent_reservoir      = reagent_res_1.rows()[0][3:7]
    Elution.reagent_reservoir    = reagent_res_1.rows()[0][9:11]
    Wash.reagent_reservoir       = reagent_res_2.rows()[0][:5]
    Ethanol.reagent_reservoir    = reagent_res_2.rows()[0][6:]
    work_destinations            = deepwell_plate.rows()[0][:Sample.num_wells]
    final_destinations           = elution_plate.rows()[0][:Sample.num_wells]

    # pipettes.
    m300 = ctx.load_instrument('p300_multi_gen2', 'left') # Load multi pipette

    #### used tip counter and set maximum tips available
    tip_track = {
        'actual': {
            tips300Fago: FIRST_TIPS_COLUMN*8,
            tips300Beads: FIRST_TIPS_COLUMN*8,
            tips300Wash: FIRST_TIPS_COLUMN*8,
            tips300Ethanol: FIRST_TIPS_COLUMN*8,
            tips300Elution: FIRST_TIPS_COLUMN*8,
            tips300Samples: FIRST_TIPS_COLUMN*8
        },
        'counts': {
            tips300Fago: 0,
            tips300Beads: 0,
            tips300Wash: 0,
            tips300Ethanol: 0,
            tips300Elution: 0,
            tips300Samples: 0
        },
        'maxes': {
            tips300Fago: 96 - (FIRST_TIPS_COLUMN * 8),
            tips300Beads: 96 - (FIRST_TIPS_COLUMN * 8),
            tips300Wash: 96 - (FIRST_TIPS_COLUMN * 8),
            tips300Ethanol: 96 - (FIRST_TIPS_COLUMN * 8),
            tips300Elution: 96 - (FIRST_TIPS_COLUMN * 8),
            tips300Samples: 96 - (FIRST_TIPS_COLUMN * 8)
        },
    }

###############################################################################

###############################################################################

    ###############################################################################
    # STEP 1 TRANSFER FAGO
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer fago
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        fago_trips = math.ceil(Fago.reagent_volume / Fago.max_volume_allowed)
        fago_volume = Fago.reagent_volume / fago_trips
        fago_transfer_vol = []
        for i in range(fago_trips):
            fago_transfer_vol.append(fago_volume + Fago.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False

        if not m300.hw_pipette['has_tip']:
            pick_up(m300, tips300Fago)
        ctx.comment('Mixing reservoir column: ' + str(Fago.col))
        # custom_mix(m300, Fago, Fago.reagent_reservoir, vol = Fago.max_volume_allowed, rounds = 3, blow_out = False, mix_height = 3, offset = 0)

        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Fago)
            for j,transfer_vol in enumerate(fago_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Fago, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from reservoir column: ' + str(Fago.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                move_vol_multi(m300, reagent = Fago, source = Fago.reagent_reservoir,
                    dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                    pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True,
                    touch_tip = False, drop_height = -38, blow_height = -38)
            ctx.comment(' ')
            ctx.comment('Mixing sample ')
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Fago.air_gap_vol_bottom) #air gap
            drop_tip(m300, tips300Fago, True)      
            
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Fago]))
    ###############################################################################
    # END STEP 1 TRANSFER FAGO
    ###############################################################################

    ###############################################################################
    # STEP 2 TRANSFER LYSIS
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer lysis
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        lysis_trips = math.ceil(Lysis.reagent_volume / Lysis.max_volume_allowed)
        lysis_volume = Lysis.reagent_volume / lysis_trips
        lysis_transfer_vol = []
        for i in range(lysis_trips):
            lysis_transfer_vol.append(lysis_volume + Lysis.disposal_volume)
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False # Original: True 
        first_mix_done = False

        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Beads)
            for j,transfer_vol in enumerate(lysis_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Lysis, multi_well_rack_area, transfer_vol * 8)
                if change_col == True or not first_mix_done: #If we switch column because there is not enough volume left in current reservoir column we mix new column
                    ctx.comment('Mixing new reservoir column: ' + str(Lysis.col))
                    custom_mix(m300, Lysis, Lysis.reagent_reservoir[Lysis.col], vol = Lysis.max_volume_allowed,
                        rounds = 10, blow_out = False, mix_height = 3, offset = 0)
                    first_mix_done = True
                else:
                    ctx.comment('Mixing reservoir column: ' + str(Lysis.col))
                    custom_mix(m300, Lysis, Lysis.reagent_reservoir[Lysis.col],
                            vol = Lysis.max_volume_allowed, rounds = 3, blow_out = False, mix_height = 3, offset = 0)
                ctx.comment('Aspirate from reservoir column: ' + str(Lysis.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                #if j!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = Lysis, source = Lysis.reagent_reservoir[Lysis.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True, touch_tip = False, drop_height = 0, blow_height=0)
            ctx.comment(' ')
            ctx.comment('Mixing sample ')
            custom_mix(m300, Lysis, location = work_destinations[i], vol =  150,
                    rounds = 20, blow_out = False, mix_height = 0, offset = 0)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Lysis.air_gap_vol_bottom) #air gap
            drop_tip(m300, tips300Beads, False)      
        tip_track['actual'][tips300Beads] = FIRST_TIPS_COLUMN * 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Beads]))
    ###############################################################################
    # END STEP 2 TRANSFER LYSIS
    ###############################################################################

    ###############################################################################
    # STEP 3 WAIT REST
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Rest for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ###############################################################################
    # END STEP 3 WAIT REST
    ###############################################################################

    ###############################################################################
    # STEP 4 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.comment(' ')
        magdeck.engage(height = mag_height)
        ctx.delay(seconds = STEPS[STEP]['wait_time'], msg = 'Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ###############################################################################
    # END STEP 4 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################

    ###############################################################################
    # STEP 5 REMOVE SUPERNATANT
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        supernatant_trips = math.ceil((Lysis.reagent_volume + VOLUME_SAMPLE) / Lysis.max_volume_allowed)
        supernatant_volume = Lysis.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Sample.disposal_volume)
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Beads)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 0.5 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                d = work_destinations[i].top(z = -5).move(Point(x = 0))
                m300.dispense(180, d, rate = 1)
                move_vol_multi(m300, reagent = Sample, source = work_destinations[i],
                        dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = True, blow_wash=True)
            d = waste.top(z = -5).move(Point(x = 0))
            m300.dispense(180, d, rate = 1)
            drop_tip(m300, tips300Beads, True)

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Beads]))
    ###############################################################################
    # END STEP 5 REMOVE SUPERNATANT
    ###############################################################################

    ###############################################################################
    # STEP 6 MAGNET OFF
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ###############################################################################
    # END STEP 6 MAGNET OFF
    ###############################################################################

    ###############################################################################
    # STEP 7 ADD WASH
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        wash_trips = math.ceil(Wash.reagent_volume / Wash.max_volume_allowed)
        wash_volume = Wash.reagent_volume / wash_trips #136.66
        wash_transfer_vol = []
        for i in range(wash_trips):
            wash_transfer_vol.append(wash_volume + Wash.disposal_volume)
        x_offset_rs = 2.5
        rinse = False # Not needed

        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Wash)
            for transfer_vol in wash_transfer_vol:
                ctx.comment('Aspirate from reservoir 1')
                [pickup_height, change_col] = calc_height(Wash, multi_well_rack_area, transfer_vol * 8)
                d = Wash.reagent_reservoir[Wash.col].top(z = -5).move(Point(x = 0))
                m300.dispense(180, d, rate = 1)
                move_vol_multi(m300, reagent = Wash, source = Wash.reagent_reservoir[Wash.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = False, drop_height=0, blow_height=0, blow_wash=True)
            d = work_destinations[i].top(z = -5).move(Point(x = 0))
            m300.dispense(180, d, rate = 1)
            custom_mix(m300, Wash, location = work_destinations[i], vol = 180,
                    rounds = 10, blow_out = False, mix_height = 3, offset = x_offset_dest - 1, mix_height_bot=-33)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Wash.air_gap_vol_bottom) #air gap
            drop_tip(m300, tips300Wash, False)
        tip_track['actual'][tips300Wash] = FIRST_TIPS_COLUMN * 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Wash]))
    ###############################################################################
    # END STEP 7 ADD WASH
    ###############################################################################

    ###############################################################################
    # STEP 8 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ####################################################################
    # END STEP 8 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################

    ###############################################################################
    # STEP 9 REMOVE SUPERNATANT
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        supernatant_trips = math.ceil(Wash.reagent_volume / Wash.max_volume_allowed)
        supernatant_volume = Wash.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Sample.disposal_volume)
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Wash)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 0.5 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                d = work_destinations[i].top(z = -5).move(Point(x = 0))
                m300.dispense(180, d, rate = 1)
                move_vol_multi(m300, reagent = Wash, source = work_destinations[i],
                    dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                    pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False, blow_wash=True)
            d = waste.top(z = -5).move(Point(x = 0))
            m300.dispense(180, d, rate = 1)
            drop_tip(m300, tips300Wash, True)

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Wash]))
    ###############################################################################
    # STEP 9 REMOVE SUPERNATANT
    ###############################################################################

    ###############################################################################
    # STEP 10 MAGNET OFF
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ###############################################################################
    # END STEP 10 MAGNET OFF
    ###############################################################################

    ###############################################################################
    # STEP 11 ADD ETHANOL
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ethanol_trips = math.ceil(Ethanol.reagent_volume / Ethanol.max_volume_allowed)
        ethanol_volume = Ethanol.reagent_volume / ethanol_trips #136.66
        ethanol_transfer_vol = []
        for i in range(ethanol_trips):
            ethanol_transfer_vol.append(ethanol_volume + Ethanol.disposal_volume)
        x_offset_rs = 1.9
        rinse = False # Not needed

        ########
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Ethanol)
            for transfer_vol in ethanol_transfer_vol:
                ctx.comment('Aspirate from Reservoir 2')
                [pickup_height, change_col] = calc_height(Ethanol, multi_well_rack_area, transfer_vol * 8)
                move_vol_multi(m300, reagent = Ethanol, source = Ethanol.reagent_reservoir[Ethanol.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 0, blow_out = True, drop_height=0, blow_height=0)
            custom_mix(m300, Ethanol, location = work_destinations[i], vol = 180,
                    rounds = 10, blow_out = False, mix_height = 3, offset = x_offset_dest - 1, mix_height_bot=-33)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Ethanol.air_gap_vol_bottom) #air gap
            drop_tip(m300, tips300Ethanol, False)
        tip_track['actual'][tips300Ethanol] = FIRST_TIPS_COLUMN * 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Ethanol]))
    ###############################################################################
    # END STEP 11 ADD ETHANOL
    ###############################################################################

    ###############################################################################
    # STEP 12 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ####################################################################
    # END STEP 12 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################

    ###############################################################################
    # STEP 13 REMOVE SUPERNATANT
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        supernatant_trips = math.ceil(Ethanol.reagent_volume / Ethanol.max_volume_allowed)
        supernatant_volume = Ethanol.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Sample.disposal_volume)
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Ethanol)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 0.5 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                d = work_destinations[i].top(z = -5).move(Point(x = 0))
                m300.dispense(180, d, rate = 1)
                move_vol_multi(m300, reagent = Sample, source = work_destinations[i],
                    dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                    pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = True, blow_wash=True)
            d = waste.top(z = -5).move(Point(x = 0))
            m300.dispense(180, d, rate = 1)
            drop_tip(m300, tips300Ethanol, True)

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Ethanol]))
    ###############################################################################
    # END STEP 13 REMOVE SUPERNATANT
    ###############################################################################

    ###############################################################################
    # STEP 14 ALLOW DRY
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')

        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Dry for ' + format(STEPS[STEP]['wait_time']) + ' seconds.') # minutes=2
        ctx.comment(' ')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)
    ###############################################################################
    # END STEP 14 ALLOW DRY
    ###############################################################################

    ###############################################################################
    # STEP 15 MAGNET OFF
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ###############################################################################
    # END STEP 15 MAGNET OFF
    ###############################################################################
    
    ###############################################################################
    # STEP 16 ADD ELUTION
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        elution_trips = math.ceil(Elution.reagent_volume / Elution.max_volume_allowed)
        elution_volume = Elution.reagent_volume / elution_trips
        elution_wash_vol = []
        for i in range(elution_trips):
            # 10 extra of elution
            elution_wash_vol.append(elution_volume + Sample.disposal_volume + 10)
        x_offset_rs = 2.5

        ########
        # Water or elution buffer
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -0.8 * find_side(i) * x_offset_rs # Original 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Elution)
            for transfer_vol in elution_wash_vol:
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Elution, multi_well_rack_area, transfer_vol*8)
                ctx.comment('Aspirate from Reservoir column: ' + str(Elution.col))
                ctx.comment('Pickup height is ' + str(pickup_height))

                move_vol_multi(m300, reagent = Elution, source = Elution.reagent_reservoir[Elution.col],
                        dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 0, blow_out = True, drop_height=0, blow_height=0)
            ctx.comment(' ')
            ctx.comment('Mixing sample with Elution')
            custom_mix(m300, Sample, work_destinations[i], vol = 40, rounds = 10,
                    blow_out = False, mix_height = 0.1, offset = x_offset_dest, mix_height_bot=-33)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            drop_tip(m300, tips300Elution, True)
        tip_track['actual'][tips300Elution] = FIRST_TIPS_COLUMN * 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Elution]))
    ###############################################################################
    # END STEP 16 ADD ELUTION
    ###############################################################################

    ###############################################################################
    # STEP 17 WAIT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ####################################################################
    # END STEP 17 WAIT
    ###############################################################################

    ###############################################################################
    # STEP 18 INCUBATE WAIT WITH MAGNET ON
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubate with magnet ON for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
    ####################################################################
    # STEP 18 WAIT FOR 5'
    ###############################################################################

    ###############################################################################
    # STEP 19 TRANSFER TO ELUTION PLATE
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        elution_trips = math.ceil(Elution.reagent_volume / Elution.max_volume_allowed)
        elution_volume = Elution.reagent_volume / elution_trips
        elution_vol = []
        for i in range(elution_trips):
            elution_vol.append(elution_volume + Elution.disposal_volume)
        x_offset_rs = 2
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300, tips300Samples)
            for transfer_vol in elution_vol:
                #Pickup_height is fixed here
                pickup_height = 0.2
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')

                move_vol_multi(m300, reagent = Elution, source = work_destinations[i],
                        dest = final_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = True, touch_tip=True)
            drop_tip(m300, tips300Samples, True)

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][tips300Samples]))

        if SET_TEMP_ON == True:
            tempdeck.set_temperature(TEMPERATURE)
    ###############################################################################
    # END STEP 19 TRANSFER TO ELUTION PLATE
    ###############################################################################

    '''if not ctx.is_simulating():
        with open(file_path,'w') as outfile:
            json.dump(STEPS, outfile)'''

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('Homing robot')
    ctx.comment('###############################################')
    ctx.comment(' ')
    ctx.home()
###############################################################################
    # Export the time log to a tsv file
    if not ctx.is_simulating():
        with open(file_path, 'w') as f:
            f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
            for key in STEPS.keys():
                row = str(key)
                for key2 in STEPS[key].keys():
                    row += '\t' + format(STEPS[key][key2])
                f.write(row + '\n')
        f.close()

    # Light flash end of program
    import os
    #os.system('mpg123 /etc/audio/speaker-test.mp3')
    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        ctx._hw_manager.hardware.set_lights(button=(0, 0 ,1))
        time.sleep(0.3)
    ctx._hw_manager.hardware.set_lights(button=(0, 1 ,0))
    ctx.comment('Finished! \nMove deepwell plate (slot 5) to Station C for MMIX addition and PCR preparation.')

    def sum_tips(tips):
        nonlocal ctx
        res = 0
        for tip in tips.keys():
            res += tips[tip]
        return res

    ctx.comment('Used tips in total: '+str(sum_tips(tip_track['counts'])))
    ctx.comment('Used racks in total: '+str(sum_tips(tip_track['counts'])/96))
    ctx.comment('Available tips: '+str(5*96))