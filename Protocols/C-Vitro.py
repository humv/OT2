import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import os
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'Station C - Vitro',
    'author': 'Aitor Gastaminza, José Luis Villanueva (Hospital Clinic Barcelona) & Alex Gasulla, Manuel Alba & Daniel Peñil',
    'source': 'Hospital Clínic Barcelona & HU Marqués de Valdecilla',
    'apiLevel': '2.3',
    'description': 'Protocol for sample setup (C) prior to qPCR (Vitro)'
    }

'''
'technician': '$technician',
'date': '$date'
'''
#Defined variables
##################
NUM_SAMPLES                 = 96 # Including controls. 94 samples + 2 controls = 96
MMIX_VOL_PER_SAMPLE         = 12
VOLUME_SAMPLE               = 8  # Volume of the sample
SET_TEMP_ON                 = True # True. Do you want to start temperature module?
TEMPERATURE                 = 4  # Temperature of temp module
##################

run_id                      = 'C_Vitro'
air_gap_vol                 = 5
air_gap_sample              = 2

# Tune variables
extra_dispensal             = 1  # Extra volume for master mix in each distribute transfer
diameter_screwcap           = 8.25  # Diameter of the screwcap
volume_cone                 = 50  # Volume in ul that fit in the screwcap cone
pipette_allowed_capacity    = 180 # Volume allowed in the pipette of 200µl
x_offset                    = [0,0]

size_transfer = math.floor(pipette_allowed_capacity / MMIX_VOL_PER_SAMPLE) # Number of wells the distribute function will fill

# Calculated variables
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on

def run(ctx: protocol_api.ProtocolContext):
    ctx.comment('Actual used columns: ' + str(num_cols))

    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': True, 'description': 'Transfer Mmix'},
        2: {'Execute': True, 'description': 'Transfer samples'},
        3: {'Execute': True, 'description': 'Transfer negative control'},
        4: {'Execute': True, 'description': 'Transfer positive control'}
    }

    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    #Folder and file_path for log time
    folder_path = '/var/lib/jupyter/notebooks' + run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/Station_C_Vitro_time_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none'):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.rinse = bool(rinse)
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.delay = delay
            self.num_wells = num_wells
            self.col = 0
            self.vol_well = 0
            self.h_cono = h_cono
            self.v_cono = v_fondo
            self.unused=[]
            self.tip_recycling = tip_recycling
            self.vol_well_original = reagent_reservoir_volume / num_wells

    # Reagents and their characteristics
    Mmix = Reagent(name = 'Mmix',
                      rinse = False,
                      flow_rate_aspirate = 3,
                      flow_rate_dispense = 3,
                      reagent_reservoir_volume = 1800,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    Samples = Reagent(name='Samples',
                      rinse=False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume=50,
                      delay=0,
                      num_wells=num_cols,  # num_cols comes from available columns
                      h_cono=0,
                      v_fondo=0
                      )

    Mmix.vol_well = Mmix.vol_well_original
    Samples.vol_well = Samples.vol_well_original

    ##################
    # Custom functions
    def divide_volume(volume,max_vol):
        num_transfers=math.ceil(volume/max_vol)
        vol_roundup=math.ceil(volume/num_transfers)
        last_vol = volume - vol_roundup*(num_transfers-1)
        vol_list = [vol_roundup for v in range(1,num_transfers)]
        vol_list.append(last_vol)
        return vol_list

    def divide_destinations(l, n):
        # Divide the list of destinations in size n lists.
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def distribute_custom(pipette, volume, src, dest, waste_pool, pickup_height, extra_dispensal, dest_x_offset, disp_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        pipette.aspirate((len(dest) * volume) +
                         extra_dispensal, src.bottom(pickup_height))
        pipette.touch_tip(speed=20, v_offset=-5)
        pipette.move_to(src.top(z=5))
        pipette.aspirate(5)  # air gap
        for d in dest:
            pipette.dispense(5, d.top())
            drop = d.top(z = disp_height).move(Point(x = dest_x_offset))
            pipette.dispense(volume, drop)
            pipette.move_to(d.top(z=5))
            pipette.aspirate(5)  # air gap
        try:
            pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
        except:
            pipette.blow_out(waste_pool.bottom(pickup_height + 3))
        return (len(dest) * volume)

    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, x_offset,
                       pickup_height, rinse, disp_height, blow_out, touch_tip):
        '''
        x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        pickup_height: height from bottom where volume
        rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        blow_out, touch_tip: if True they will be done after dispensing
        '''
        # Rinse before aspirating
        if rinse == True:
            custom_mix(pipet, reagent, location = source, vol = vol,
                       rounds = 2, blow_out = True, mix_height = 0,
                       x_offset = x_offset)
        # SOURCE
        s = source.bottom(pickup_height).move(Point(x = x_offset[0]))
        pipet.aspirate(vol, s)  # aspirate liquid
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pipet.aspirate(air_gap_vol, source.top(z = -2),
                           rate = reagent.flow_rate_aspirate)  # air gap
        # GO TO DESTINATION
        drop = dest.top(z = disp_height).move(Point(x = x_offset[1]))
        pipet.dispense(vol + air_gap_vol, drop,
                       rate = reagent.flow_rate_dispense)  # dispense all
        ctx.delay(seconds = reagent.delay) # pause for x seconds depending on reagent
        if blow_out == True:
            pipet.blow_out(dest.top(z = -2))
        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -5, radius=0.5)


    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 3):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        if mix_height == 0:
            mix_height = 3
        pipet.aspirate(1, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
        for _ in range(rounds):
            pipet.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
            pipet.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        pipet.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        if blow_out == True:
            pipet.blow_out(location.top(z=-2))  # Blow out

    def calc_height(reagent, cross_section_area, aspirate_volume, min_height=0.5):
        nonlocal ctx
        ctx.comment('Remaining volume ' + str(reagent.vol_well) +
                    '< needed volume ' + str(aspirate_volume) + '?')
        if reagent.vol_well < aspirate_volume:
            reagent.unused.append(reagent.vol_well)
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
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Calculated height is ' + str(height))
            if height < min_height:
                height = min_height
            ctx.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    ####################################
    # load labware and modules
    # 24 well rack
    tuberack = ctx.load_labware(
        'opentrons_24_aluminumblock_generic_2ml_screwcap', '2',
        'Opentrons 24 Well Aluminum Block with Generic 2 mL Screwcap')

    ############################################
    # tempdecks
    tempdeck_orig = ctx.load_module('Temperature Module Gen2', '4')
    tempdeck_dest = ctx.load_module('Temperature Module Gen2', '1')

    if(SET_TEMP_ON):
        tempdeck_orig.set_temperature(TEMPERATURE)
        tempdeck_dest.set_temperature(TEMPERATURE)

    ##################################
    # Sample plate - comes from B
    source_plate = tempdeck_orig.load_labware(
        'kingfisher_96_aluminumblock_200ul', 
        'Kingfisher 96 Aluminum Block 200 uL')
    samples = source_plate.wells()[:NUM_SAMPLES]

    ##################################
    # qPCR plate - final plate, goes to PCR
    qpcr_plate = tempdeck_dest.load_labware(
        'opentrons_96_aluminumblock_nest_wellplate_100ul', 
        'Opentrons 96 Well Aluminum Block with NEST Well Plate 100 uL')

    ##################################
    # Load Tipracks
    tips20 = [
        ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
        for slot in ['5']
    ]

    tips200 = [
        ctx.load_labware('opentrons_96_filtertiprack_200ul', slot)
        for slot in ['3']
    ]

    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    Mmix.reagent_reservoir = tuberack.rows()[0][0] # A1

    # setup up sample sources and destinations
    samples = source_plate.wells()[2:NUM_SAMPLES]
    pcr_wells = qpcr_plate.wells()[:NUM_SAMPLES]
    pcr_wells_samples = qpcr_plate.wells()[2:NUM_SAMPLES]

    # Divide destination wells in small groups for P300 pipette
    dests = list(divide_destinations(pcr_wells, size_transfer))

    # pipettes
    p20 = ctx.load_instrument(
        'p20_single_gen2', mount='right', tip_racks=tips20)
    p300 = ctx.load_instrument(
        'p300_single_gen2', mount='left', tip_racks=tips200)

    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {p300: 0,
                   p20: 0},
        'maxes': {p300: 96 * len(p300.tip_racks),
                    p20: 96 * len(p20.tip_racks)}
    }

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        if not ctx.is_simulating():
            if tip_track['counts'][pip] == tip_track['maxes'][pip]:
                ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                pip.reset_tipracks()
                tip_track['counts'][pip] = 0

        if not pip.hw_pipette['has_tip']:
            pip.pick_up_tip()
    ##########

    ############################################################################
    # STEP 1: TRANSFER MMIX
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        pick_up(p300)
        used_vol = []

        for dest in dests:
            aspirate_volume = MMIX_VOL_PER_SAMPLE * len(dest) + extra_dispensal
            used_vol_temp = distribute_custom(p300, volume = MMIX_VOL_PER_SAMPLE,
                src = Mmix.reagent_reservoir, dest = dest,
                waste_pool = Mmix.reagent_reservoir, pickup_height = 0.2,
                extra_dispensal = extra_dispensal, dest_x_offset = 2, disp_height = -1)
            used_vol.append(used_vol_temp)
        p300.drop_tip(home_after = False)
        tip_track['counts'][p300]+=1

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 2: TRANSFER SAMPLES
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        for s, d in zip(samples, pcr_wells_samples):
            pick_up(p20)
            move_vol_multichannel(p20, reagent = Samples, source = s, dest = d,
                    vol = VOLUME_SAMPLE, air_gap_vol = air_gap_sample, x_offset = x_offset,
                    pickup_height = 0.2, disp_height = 0, rinse = False,
                    blow_out=True, touch_tip=True)
            p20.drop_tip(home_after = False)
            tip_track['counts'][p20]+=1

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3: TRANSFER NEGATIVE CONTROL
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        pick_up(p20)
        s = tuberack.rows()[0][1]   # A2
        d = qpcr_plate.wells()[1]   # B1
        move_vol_multichannel(p20, reagent = Samples, source = s, dest = d,
                vol = VOLUME_SAMPLE, air_gap_vol = air_gap_sample, x_offset = x_offset,
                pickup_height = 0.2, disp_height = 0, rinse = False,
                blow_out=True, touch_tip=True)
        p20.drop_tip(home_after = False)
        tip_track['counts'][p20]+=1

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 4: TRANSFER POSITIVE CONTROL
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        pick_up(p20)
        s = tuberack.rows()[0][2]   # A3
        d = qpcr_plate.wells()[0]   # A1
        move_vol_multichannel(p20, reagent = Samples, source = s, dest = d,
                vol = VOLUME_SAMPLE, air_gap_vol = air_gap_sample, x_offset = x_offset,
                pickup_height = 0.2, disp_height = 0, rinse = False,
                blow_out=True, touch_tip=True)
        p20.drop_tip(home_after = False)
        tip_track['counts'][p20]+=1

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


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

    ############################################################################
    # Light flash end of program
    for i in range(3):
        #ctx._hw_manager.hardware.set_lights(rails=False)
        ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))
        time.sleep(0.3)
        #ctx._hw_manager.hardware.set_lights(rails=True)
        ctx._hw_manager.hardware.set_lights(button=(0, 0 ,1))
        time.sleep(0.3)
    ctx._hw_manager.hardware.set_lights(button=(0, 1 ,0))
    ctx.comment('Finished! \nMove plate to PCR')

    total_used_vol = np.sum(used_vol)
    total_needed_volume = total_used_vol
    ctx.comment('Total Mmix used volume is: ' + str(total_used_vol) + '\u03BCl.')
    ctx.comment('Needed Mmix volume is ' +
                str(total_needed_volume + extra_dispensal*len(dests)) +'\u03BCl')
    ctx.comment('Mmix remaining in tubes is: ' +
                format(np.sum(Mmix.unused) + extra_dispensal * len(dests) + Mmix.vol_well) + '\u03BCl.')
    ctx.comment('200 ul Used tips in total: ' + str(tip_track['counts'][p300]))
    ctx.comment('200 ul Used racks in total: ' + str(tip_track['counts'][p300] / 96))
    ctx.comment('20 ul Used tips in total: ' + str(tip_track['counts'][p20]))
    ctx.comment('20 ul Used racks in total: ' + str(tip_track['counts'][p20] / 96))
