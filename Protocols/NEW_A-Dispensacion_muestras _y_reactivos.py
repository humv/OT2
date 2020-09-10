import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import os
from timeit import default_timer as timer
from datetime import datetime
import subprocess

# metadata
metadata = {
    'protocolName': 'Station A - Sample dispensing',
    'author': 'Aitor Gastaminza, Alex Gasulla & José Luis Villanueva (Hospital Clinic Barcelona),  Manuel Alba & Daniel Peñil',
    'source': 'Hospital Clínic Barcelona & HU Marqués de Valdecilla',
    'apiLevel': '2.6',
    'description': 'Protocol for sample dispensing'
}

'''
'technician': '$technician',
'date': '$date'
'''
 
################################################
# CHANGE THESE VARIABLES ONLY
################################################
NUM_CONTROL_SPACES      = 2  # The control spaces are being ignored at the first cycles
NUM_REAL_SAMPLES        = 94   
NUM_MIXES               = 0
VOLUME_SAMPLE           = 200 # Sample volume to place in deepwell

SOUND_NUM_PLAYS         = 3
PHOTOSENSITIVE          = False # True if it has photosensitive reagents
LYSIS_VOLUME_PER_SAMPLE = 265 # ul per sample.
BEADS_VOLUME_PER_SAMPLE = 10 # ul per sample.
PK_VOLUME_PER_SAMPLE    = 10 # ul per sample.
################################################

recycle_tip             = True
num_samples             = NUM_CONTROL_SPACES + NUM_REAL_SAMPLES
air_gap_vol_sample      = 25
extra_dispensal         = 1
run_id                  = 'preparacion_tipo_A'
path_sounds             = '/var/lib/jupyter/notebooks/sonidos/'
sonido_defecto          = 'finalizado.mp3'
volume_mix              = 500 # Volume used on mix
x_offset                = [0,0]

lysys_pipette_capacity  = 900 # Volume allowed in the pipette of 1000µl
size_transfer           = math.floor(lysys_pipette_capacity / LYSIS_VOLUME_PER_SAMPLE) # Number of wells the distribute function will fill


def run(ctx: protocol_api.ProtocolContext):
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description and times
        1: {'Execute': True, 'description': 'Dispensar Lysys'},
        2: {'Execute': True, 'description': 'Mezclar y dispensar muestras ('+str(VOLUME_SAMPLE)+'ul)'}
    }
    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    #Folder and file_path for log time
    if not ctx.is_simulating():
        folder_path = '/var/lib/jupyter/notebooks/' + run_id
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/StationA_time_log.txt'

    # Define Reagents as objects with their properties
    class Simple_Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, delay):
            self.name               = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.delay              = delay 

    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense,  
            flow_rate_aspirate_mix, flow_rate_dispense_mix, air_gap_vol_bottom, air_gap_vol_top, disposal_volume, max_volume_allowed, reagent_volume, placed_in_multi, v_fondo):
            self.name               = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = 0.5
            self.flow_rate_dispense_mix = 0.5
            self.air_gap_vol_bottom = 5
            self.air_gap_vol_top = 0
            self.disposal_volume = 1
            self.max_volume_allowed = 18
            self.reagent_volume = BEADS_VOLUME_PER_SAMPLE
            self.placed_in_multi = True
            self.v_fondo = 695 #1.95 * multi_well_rack_area / 2, #Prismatic

    # Reagents and their characteristics
    Samples = Simple_Reagent(name                  = 'Samples',
                      flow_rate_aspirate    = 25,
                      flow_rate_dispense    = 100,
                      delay                 = 0
                      ) 

    Beads = Reagent(name = 'Beads',
                    flow_rate_aspirate = 0.5,
                    flow_rate_dispense = 0.5,
                    flow_rate_aspirate_mix = 0.5,
                    flow_rate_dispense_mix = 0.5,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    max_volume_allowed = 18,
                    reagent_volume = BEADS_VOLUME_PER_SAMPLE,
                    placed_in_multi = True,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    Lysis = Simple_Reagent(name                      = 'Lysis',
                     flow_rate_aspirate        = 50,
                     flow_rate_dispense        = 100,
                     delay                     = 0
                     ) 

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('CONTROLES: ' + str(NUM_CONTROL_SPACES))  
    ctx.comment('MUESTRAS: ' + str(NUM_REAL_SAMPLES)) 
    ctx.comment('###############################################')
    ctx.comment(' ')

    ##################
    # Custom functions
    ##################

    
    def log_step_start():
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('PASO '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        return datetime.now()

    def log_step_end(start):
        end = datetime.now()
        time_taken = (end - start)
        STEPS[STEP]['Time:'] = str(time_taken)

        ctx.comment(' ')
        ctx.comment('Paso ' + str(STEP) + ': ' +STEPS[STEP]['description'] + ' hizo un tiempo de ' + str(time_taken))
        ctx.comment(' ')

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
            #pipet.blow_out(dest.top(z = -2))
            pipet.blow_out(dest.top(z = disp_height))
        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -10)

        if air_gap_vol != 0:
            # pipet.move_to(dest.top(z = disp_height))
            pipet.air_gap(air_gap_vol, height = disp_height) #air gap
    
    def distribute_custom(pipette, reagent, volume, src, dest, waste_pool, pickup_height, extra_dispensal, dest_x_offset, disp_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        pipette.aspirate((len(dest) * volume) +extra_dispensal
                         , src.bottom(pickup_height), rate = reagent.flow_rate_aspirate)
        pipette.move_to(src.top(z=5))
        pipette.aspirate(air_gap_vol_sample, rate = reagent.flow_rate_aspirate)  # air gap
        for d in dest:
            pipette.dispense(volume + air_gap_vol_sample, d.top(), rate = reagent.flow_rate_dispense)
            pipette.move_to(d.top(z=5))
            pipette.aspirate(air_gap_vol_sample, rate = reagent.flow_rate_dispense)  # air gap
        try:
            pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
        except:
            pipette.blow_out(waste_pool.top(pickup_height + 3))
        return (len(dest) * volume)

    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 5):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        if mix_height <= 0:
            mix_height = 3

        pipet.aspirate(1, location = location.bottom(
                        z = source_height).move(Point(x = x_offset[0])), rate = reagent.flow_rate_aspirate)

        for _ in range(rounds):
            pipet.aspirate(vol, location = location.bottom(
                z = source_height).move(Point(x = x_offset[0])), rate = reagent.flow_rate_aspirate)
            pipet.dispense(vol, location = location.bottom(
                z = mix_height).move(Point(x = x_offset[1])), rate = reagent.flow_rate_dispense)

        pipet.dispense(1, location = location.bottom(
            z = mix_height).move(Point(x = x_offset[1])), rate = reagent.flow_rate_dispense)

        if blow_out == True:
            pipet.blow_out(location.top(z = -2))  # Blow out

    def generate_source_table(source):
        '''
        Concatenate the wells frome the different origin racks
        '''
        num_cols = math.ceil(num_samples / 8)
        s = []
        for i  in range(num_cols):
            if i < 6:
                s += source[0].columns()[i] + source[1].columns()[i]
            else:
                s += source[2].columns()[i - 6] + source[3].columns()[i - 6]
        return s

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
        pip.pick_up_tip()
    
    def pick_up_tip(pip, tips):
        nonlocal tip_track
        #if not ctx.is_simulating():
        if recycle_tip:
            pip.pick_up_tip(tips[0].wells()[0])
        else:
            if tip_track['counts'][pip] >= tip_track['maxes'][pip]:
                for i in range(3):
                    ctx._hw_manager.hardware.set_lights(rails=False)
                    ctx._hw_manager.hardware.set_lights(button=(1, 0 ,0))
                    time.sleep(0.3)
                    ctx._hw_manager.hardware.set_lights(rails=True)
                    ctx._hw_manager.hardware.set_lights(button=(0, 0 ,1))
                    time.sleep(0.3)
                ctx._hw_manager.hardware.set_lights(button=(0, 1 ,0))
                ctx.pause('Cambiar ' + str(pip.max_volume) + ' µl tipracks antes del pulsar Resume.')
                pip.reset_tipracks()
                tip_track['counts'][pip] = 0
                tip_track['num_refills'][pip] += 1
            pip.pick_up_tip()

    def run_quiet_process(command):
        subprocess.check_output('{} &> /dev/null'.format(command), shell=True)

    def play_sound(filename):
        print('Speaker')
        print('Next\t--> CTRL-C')
        try:
            run_quiet_process('mpg123 {}'.format(path_sounds + filename + '.mp3'))
            run_quiet_process('mpg123 {}'.format(path_sounds + sonido_defecto))
            run_quiet_process('mpg123 {}'.format(path_sounds + filename + '.mp3'))

        except KeyboardInterrupt:
            pass
            print()
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

    def finish_run():
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
        ctx._hw_manager.hardware.set_lights(button = True, rails =  False)

        used_tips = tip_track['num_refills'][p1000] * 96 * len(p1000.tip_racks) + tip_track['counts'][p1000]
        ctx.comment('Puntas de 1000 ul utilizadas: ' + str(used_tips) + ' (' + str(round(used_tips / 96, 2)) + ' caja(s))')
        ctx.comment('###############################################')

        if not ctx.is_simulating():
            for i in range(SOUND_NUM_PLAYS):
                if i > 0:
                    time.sleep(60)
                play_sound('finished_process_esp')

            return finish_time
    
    def divide_destinations(l, n):
        # Divide the list of destinations in size n lists.
        for i in range(0, len(l), n):
            yield l[i:i + n]

    ####################################
    # load labware and modules

    ####################################
    # Load Sample racks
    if num_samples <= 48:
        rack_num = 2
        ctx.comment('Used source racks are ' + str(rack_num))
    else:
        rack_num = 4

    source_racks = [ctx.load_labware(
        'opentrons_24_tuberack_nest_2ml_snapcap', slot,
        'source tuberack with snapcap' + str(i + 1)) for i, slot in enumerate(['4', '1', '5', '2'][:rack_num])
    ]

    lysys_rack = ctx.load_labware('opentrons_6_tuberack_falcon_50ml_conical', '3','source tuberack with snapcap')
    reagent_res = ctx.load_labware('nest_12_reservoir_15ml', '8','reagent deepwell plate')

    ##################################
    # Destination plate
    dest_plate = ctx.load_labware(
        'nest_96_wellplate_2ml_deep', '6',
        'NEST 96 Deepwell Plate 2mL')

    ####################################
    # Load tip_racks
    tips1000 = [ctx.load_labware(
        'opentrons_96_filtertiprack_1000ul', slot, 
        '1000µl filter tiprack') for slot in ['7']]

    tips20 = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot, ' filter tiprack')
        for slot in ['11']]

    ################################################################################
    # setup samples and destinations
    sample_sources_full = generate_source_table(source_racks)
    sample_sources      = sample_sources_full[NUM_CONTROL_SPACES:num_samples]
    destinations        = dest_plate.wells()[NUM_CONTROL_SPACES:num_samples]
    lysys_source        = lysys_rack.wells_by_name()['B3']
    dests_lysis         = list(divide_destinations(destinations, size_transfer))

    p1000 = ctx.load_instrument(
        'p1000_single_gen2', 'right', 
        tip_racks = tips1000) # load P1000 pipette

    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {p1000: 0},
        'maxes': {p1000: 96 * len(p1000.tip_racks)}, #96 tips per tiprack * number or tipracks in the layout
        'num_refills' : {p1000 : 0},
        'tips': { p1000: [tip for rack in tips1000 for tip in rack.rows()[0]]}

    }


    start_run()
    ############################################################################
    # STEP 1: ADD LYSIS 
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = log_step_start()

        used_vol = []
        for dest in dests_lysis:
            if not p1000.hw_pipette['has_tip']:
                 pick_up_tip(p1000,tips1000)
            used_vol_temp = distribute_custom(p1000, Lysis, volume = LYSIS_VOLUME_PER_SAMPLE,
                src = lysys_source, dest = dest,
                waste_pool = lysys_source, pickup_height = 1,
                extra_dispensal = extra_dispensal, dest_x_offset = 2, disp_height = -1)
            used_vol.append(used_vol_temp)

        p1000.drop_tip(home_after = False)

        log_step_end(start)

    ############################################################################
    # STEP 2: MIX AND MOVE SAMPLES
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        start = datetime.now()
        for s, d in zip(sample_sources, destinations):
            if not p1000.hw_pipette['has_tip']:
                pick_up(p1000)

            # Mix the sample BEFORE dispensing
            if NUM_MIXES > 0:
                custom_mix(p1000, reagent = Samples, location = s, vol = volume_mix, 
                    rounds = NUM_MIXES, blow_out = True, mix_height = 15, x_offset = x_offset)

            move_vol_multichannel(p1000, reagent = Samples, source = s, dest = d,
                vol = VOLUME_SAMPLE, air_gap_vol = air_gap_vol_sample, x_offset = x_offset,
                pickup_height = 3, rinse = False, disp_height = -10,
                blow_out = True, touch_tip = False)

            p1000.drop_tip(home_after = False)
            tip_track['counts'][p1000] += 1

        # Time statistics
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] +
                    ' took ' + str(time_taken))
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
    # from opentrons.drivers.rpi_drivers import gpio

    finish_run()