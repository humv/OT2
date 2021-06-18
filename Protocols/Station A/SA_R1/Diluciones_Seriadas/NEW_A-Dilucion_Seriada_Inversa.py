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
    'protocolName': 'Station A - Sample dilucion',
    'author': 'CIC',
    'source': 'HU Marqués de Valdecilla',
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
NUM_REAL_SAMPLES        = 88
NUM_MIX_SPACES          = 8     # The control spaces are being ignored at the last cycles

VOLUME_ANTBIOTIC        = 100   # Sample volume to place in deepwell
VOLUME_SAMPLE           = 50   # Sample volume to place in deepwell

BEADS_WELL_FIRST_TIME_NUM_MIXES     = 20
BEADS_WELL_NUM_MIXES            = 3

TUBE_NUM_MIXES          = 0

SOUND_NUM_PLAYS         = 1
PHOTOSENSITIVE          = False # True if it has photosensitive reagents

DEFAULT_DEAD_VOL            = 700

################################################

air_gap_vol_sample      = 25
run_id                  = 'A-Dilucion_Seriada_Inversa'
path_sounds             = '/var/lib/jupyter/notebooks/sonidos/'
sonido_defecto          = 'finalizado.mp3'
x_offset                = [0,0]
switch_off_lights       = False # Switch of the lights when the program finishes
multi_well_rack_area        = 8 * 71 #Cross section of the 12 well reservoir

num_cols_antibiotic = 1
num_cols_caldo = 11

def run(ctx: protocol_api.ProtocolContext):
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description and times
        1: {'Execute': True, 'description': 'Mezclar cultivos'}
        
    }
    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    #Folder and file_path for log time
    if not ctx.is_simulating():
        folder_path = '/var/lib/jupyter/notebooks/' + run_id
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/time_log.txt'

    # Define Reagents as objects with their properties
        #Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, flow_rate_aspirate_mix, flow_rate_dispense_mix,
        air_gap_vol_bottom, air_gap_vol_top, disposal_volume, rinse, max_volume_allowed, reagent_volume, reagent_reservoir_volume, num_wells, h_cono, v_fondo, tip_recycling = 'none', dead_vol = DEFAULT_DEAD_VOL):
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

    # Reagents and their characteristics
    Samples = Reagent(name = 'Antibioticos',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 25,
                    flow_rate_dispense_mix = 50,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = VOLUME_SAMPLE,
                    reagent_reservoir_volume = NUM_REAL_SAMPLES * VOLUME_SAMPLE * 1.1,
                    num_wells = math.ceil(NUM_REAL_SAMPLES  * VOLUME_SAMPLE * 1.1 / 11500),
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    
    Caldo = Reagent(name = 'Caldo ',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 25,
                    flow_rate_dispense_mix = 50,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = VOLUME_SAMPLE,
                    reagent_reservoir_volume = NUM_REAL_SAMPLES * VOLUME_SAMPLE * 1.1,
                    num_wells = math.ceil(NUM_REAL_SAMPLES  * VOLUME_SAMPLE * 1.1 / 11500),
                    h_cono = 1.95,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('VALORES DE VARIABLES')
    ctx.comment(' ')
    ctx.comment('Número de muestras: ' + str(NUM_REAL_SAMPLES) + ' (' + str(num_cols_antibiotic) + ' columnas)')
    ctx.comment(' ')
    ctx.comment('Volumen de muestra a mover al deepwell: ' + str(VOLUME_SAMPLE) + ' ul')
    ctx.comment(' ')
    ctx.comment('Número de mezclas en la muestra: ' + str(TUBE_NUM_MIXES))
    ctx.comment(' ')
    ctx.comment('Repeticiones del sonido final: ' + str(SOUND_NUM_PLAYS))
    ctx.comment('Foto-sensible: ' + str(PHOTOSENSITIVE))
    ctx.comment(' ')

    ##################
    # Custom functions
    def move_multichanel_diluido(dest):
        beads_trips = math.ceil(Caldo.reagent_volume / Caldo.max_volume_allowed)
        beads_volume = Caldo.reagent_volume / beads_trips #136.66
        beads_transfer_vol = []
        for i in range(beads_trips):
            beads_transfer_vol.append(beads_volume + Caldo.disposal_volume)

        x_offset_source = 0
        x_offset_dest   = 0
        rinse = False # Original: True
        first_mix_done = False
        
        if not m300.hw_pipette['has_tip']:
            pick_up_multi(m300)
        for i in range(num_cols_caldo):
            ctx.comment("Column: " + str(i))
            for j,transfer_vol in enumerate(beads_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Caldo, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from reservoir column: ' + str(Caldo.col))
                ctx.comment('Pickup height is ' + str(round(pickup_height, 2)) + ' mm')
 
                move_vol_multi(m300, reagent = Caldo, source = Caldo.reagent_reservoir,
                        dest = dest[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                        pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 2, blow_out = True, touch_tip = True, drop_height = -1)
        m300.drop_tip(home_after = False)
        tip_track['counts'][m300] += 8
        
    
        

        
    
    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, x_offset,
                       pickup_height, rinse, disp_height, blow_out, touch_tip):
        '''
        x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        pickup_height: height from bottom where volume
        rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        blow_out, touch_tip: if True they will be done after dispensing
        '''

        # SOURCE
        s = source.bottom(pickup_height).move(Point(x = x_offset[0]))
        pipet.aspirate(vol, s, rate = reagent.flow_rate_aspirate)  # aspirate liquid
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pipet.aspirate(air_gap_vol, source.top(z = -2),
                           rate = reagent.flow_rate_aspirate)  # air gap

        # GO TO DESTINATION
        drop = dest.top(z = disp_height).move(Point(x = x_offset[1]))
        pipet.dispense(vol + air_gap_vol, drop,
                       rate = reagent.flow_rate_dispense)  # dispense all

        #ctx.delay(seconds = reagent.delay) # pause for x seconds depending on reagent

        if blow_out == True:
            pipet.blow_out(dest.top(z = disp_height))

        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -10)

        if air_gap_vol != 0:
            pipet.air_gap(air_gap_vol, height = disp_height) #air gap

    
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

    def pick_up_multi(pip):
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

        used_tips_1000 = tip_track['num_refills'][p1000] * 96 * len(p1000.tip_racks) + tip_track['counts'][p1000]
        used_tips_300 = tip_track['num_refills'][m300] * 96 * len(m300.tip_racks) + tip_track['counts'][m300]
        ctx.comment('Puntas de 1000 ul utilizadas: ' + str(used_tips_1000) + ' (' + str(round(used_tips_1000 / 96, 2)) + ' caja(s))')
        ctx.comment('Puntas de 300 ul utilizadas: ' + str(used_tips_300) + ' (' + str(round(used_tips_300 / 96, 2)) + ' caja(s))')
        ctx.comment('###############################################')

        if not ctx.is_simulating():
            for i in range(SOUND_NUM_PLAYS):
                if i > 0:
                    time.sleep(60)
                play_sound('finished_process_esp')

            return finish_time

    def validate_parameters():
        result = True

        if NUM_REAL_SAMPLES + NUM_REAL_SAMPLES > 96:
            ctx.comment ("ERROR: La suma de cada gradilla del número de muestras ("+ str(NUM_REAL_SAMPLES) +")  es mayor que 96, verifique los valores y vuelva a cargar el protocolo.")
            result = False

        return result

    ####################################
    # load labware and modules

    ####################################
    # Load Sample racks
    
    source_antibiotic =  ctx.load_labware(
        'opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap', '2',
        'Opentrons 24 Tuberack Eppendorf 2ml Safelock Snapcap')

    source_caldo = ctx.load_labware(
        'nest_1_reservoir_195ml', '1',
        'NEST 1 Reservoir 195ul')

    source_diluido = ctx.load_labware(
        'nest_1_reservoir_195ml', '4',
        'NEST 1 Reservoir diluido 195ul')

    ##################################
    # Destination plate
    # Destination
    dest_plate1 = ctx.load_labware(
        'nest_96_wellplate_200ul_flat', '3',
        'NEST 96 Well Plate 200ul 1')
        
    dest_plate2 = ctx.load_labware(
        'nest_96_wellplate_200ul_flat', '5',
        'NEST 96 Well Plate 200ul 2')
        
    dest_plate3 = ctx.load_labware(
        'nest_96_wellplate_200ul_flat', '7',
        'NEST 96 Well Plate 200ul 3')


    

    ####################################
    # Load tip_racks

    tips200 = [
        ctx.load_labware('opentrons_96_filtertiprack_200ul', slot)
        for slot in ['10']
    ]

    tips1000 = [
        ctx.load_labware('opentrons_96_filtertiprack_1000ul', slot)
        for slot in ['11']
    ]
    ################################################################################
    # setup samples and destinations
    Caldo.reagent_reservoir      = source_diluido.wells()[0]
    destinations1        = dest_plate1.rows()[0][1:]
    destinations1.reverse()
    destinations2        = dest_plate2.rows()[0][1:]
    destinations2.reverse()
    destinations3        = dest_plate3.rows()[0][1:]
    destinations3.reverse()

    m300 = ctx.load_instrument(
        'p300_multi_gen2', 'right', 
        tip_racks = tips200) # load P1000 pipette

    p1000 = ctx.load_instrument(
        'p1000_single_gen2', mount = 'left', tip_racks = tips1000)
    
    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {m300: 0},
        'maxes': {m300: 96 * len(m300.tip_racks)},
        'num_refills' : {m300: 0},
        'tips': { m300: [tip for rack in tips200 for tip in rack.rows()[0]]}
    }

    if validate_parameters():

        start_run()

        ############################################################################
        # STEP 1: PREPARE DISPENSE AND MOVE SAMPLES
        ############################################################################
        STEP += 1
        if STEPS[STEP]['Execute'] == True:
            ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
            ctx.comment('###############################################')

            start = datetime.now()
            if not m300.hw_pipette['has_tip']:
                pick_up_multi(m300)
            move_multichanel_diluido(destinations1)
            move_multichanel_diluido(destinations2)
            move_multichanel_diluido(destinations3)           


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

        finish_run(switch_off_lights)