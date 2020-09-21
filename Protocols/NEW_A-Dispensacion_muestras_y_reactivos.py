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
    'author': 'Aitor Gastaminza, Alex Gasulla & José Luis Villanueva (Hospital Clinic Barcelona),  Manuel Alba ,Daniel Peñil & David Martínez',
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
NUM_CONTROL_SPACES          = 0  # The control spaces are being ignored at the first cycles
NUM_REAL_SAMPLES            = 16

VOLUME_SAMPLE               = 200 # Sample volume to place in deepwell
LYSIS_VOLUME_PER_SAMPLE     = 265 # ul per sample.
PK_VOLUME_PER_SAMPLE        = 13 # ul per sample.
BEADS_VOLUME_PER_SAMPLE     = 13 # ul per sample.

NUM_BEFORE_MIXES            = 0
NUM_AFTER_MIXES             = 1

MAX_LYSYS_DISPENSE_PER_TIP  = 48 # max number of samples dispensed with the same lysys tip. ex: 48, means two tips used to dispense lysys to 96 samples.

SOUND_NUM_PLAYS             = 0
PHOTOSENSITIVE              = False # True if it has photosensitive reagents
################################################

recycle_tip             = False
num_samples             = NUM_CONTROL_SPACES + NUM_REAL_SAMPLES
num_cols                = math.ceil(num_samples / 8) # Columns we are working on


extra_dispensal         = 1
run_id                  = 'A-Dispensacion_muestras_y_reactivos-Magmax'
path_sounds             = '/var/lib/jupyter/notebooks/sonidos/'
sonido_defecto          = 'finalizado.mp3'
volume_mix_tuberack     = 500
volume_mix_deepwell     = (LYSIS_VOLUME_PER_SAMPLE + VOLUME_SAMPLE) * 0.75 # Volume used on mix
x_offset                = [0,0]
switch_off_lights       = False # Switch of the lights when the program finishes

lysys_pipette_capacity  = 900 # Volume allowed in the pipette of 1000µl
size_transfer           = math.floor(lysys_pipette_capacity / LYSIS_VOLUME_PER_SAMPLE) # Number of wells the distribute function will fill
multi_well_rack_area    = 8 * 71    #Cross section of the 12 well reservoir
next_well_index         = 0         # First reagent well to use

def run(ctx: protocol_api.ProtocolContext):
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description and times
        1: {'Execute': True, 'description': 'Dispensar Lysys'},
        2: {'Execute': True, 'description': 'Mezclar y dispensar muestras ('+str(VOLUME_SAMPLE)+'ul)'},
        3: {'Execute': True, 'description': 'Transferir proteinasa K ('+str(PK_VOLUME_PER_SAMPLE)+'ul)'},
        4: {'Execute': True, 'description': 'Transferir bolas magnéticas ('+str(BEADS_VOLUME_PER_SAMPLE)+'ul)'}
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
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense,flow_rate_aspirate_mix, flow_rate_dispense_mix, delay, air_gap_vol_bottom = 2, air_gap_vol_top = 0):
            self.name               = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
            self.flow_rate_dispense_mix = flow_rate_dispense_mix
            self.delay              = delay 
            self.air_gap_vol_bottom = air_gap_vol_bottom
            self.air_gap_vol_top = air_gap_vol_top
    
    def str_rounded(num):
        return str(int(num + 0.5))

    class Reagent:
        def calc_vol_well(self):
            global num_cols

            if(self.name == 'Sample'):
                self.num_wells = num_cols
                return VOLUME_SAMPLE
            elif self.placed_in_multi:
                trips = math.ceil(self.reagent_volume / self.max_volume_allowed)
                vol_trip = self.reagent_volume / trips * 8
                max_trips_well = math.floor(18000 / vol_trip)
                total_trips = num_cols * trips
                self.num_wells = math.ceil(total_trips / max_trips_well)
                return math.ceil(total_trips / self.num_wells) * vol_trip + self.dead_vol
            else:
                self.num_wells = 1
                return self.reagent_volume * num_samples

        def set_first_well(self, first_well_pos = None):
            global next_well_index
            if first_well_pos is not None and first_well_pos > next_well_index:
                self.first_well = first_well_pos
            else:
                self.first_well = next_well_index + 1

            next_well_index = self.first_well - 1 + self.num_wells
            
            return self.first_well

        def comment_vol_info(self):
            ctx.comment(self.name + ': ' + str(self.num_wells) +  (' canal' if self.num_wells == 1 else ' canales') + ' desde el canal '+ str(self.first_well) +' en el reservorio de 12 canales con un volumen de ' + str_rounded(self.vol_well_original) + ' uL cada uno')
        
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense,  
                air_gap_vol_bottom, disposal_volume, max_volume_allowed, reagent_volume, v_fondo, 
                flow_rate_aspirate_mix = 0.5, flow_rate_dispense_mix = 0.5, air_gap_vol_top = 0, 
                dead_vol = 700, first_well = None, placed_in_multi = False):
            self.name               = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
            self.flow_rate_dispense_mix = flow_rate_aspirate_mix
            self.air_gap_vol_bottom = 2
            self.air_gap_vol_top = air_gap_vol_top
            self.disposal_volume = 1
            self.max_volume_allowed = 18
            self.col = 0
            self.reagent_volume = BEADS_VOLUME_PER_SAMPLE
            self.v_fondo = 695 #1.95 * multi_well_rack_area / 2, #Prismatic
            self.v_cono = v_fondo
            self.dead_vol = dead_vol
            self.placed_in_multi = placed_in_multi
            self.vol_well_original = self.calc_vol_well() if reagent_volume * num_samples > 0 else 0
            self.first_well = self.set_first_well(first_well)
            self.vol_well = self.vol_well_original
            self.delay = 0

    # Reagents and their characteristics
    Samples = Simple_Reagent(name = 'Samples',
                      flow_rate_aspirate    = 25,
                      flow_rate_dispense    = 100,
                      flow_rate_aspirate_mix = 0.5,
                      flow_rate_dispense_mix = 0.5,
                      delay                 = 0,
                      air_gap_vol_bottom    = 25
                      ) 
                
    Pk = Reagent(name = 'Pk',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 0.5,
                    flow_rate_dispense_mix = 0.5,
                    air_gap_vol_bottom = 1,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    max_volume_allowed = 18,
                    reagent_volume = PK_VOLUME_PER_SAMPLE,
                    placed_in_multi = True,
                    first_well = 1,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic

    Beads = Reagent(name = 'Beads',
                    flow_rate_aspirate = 25,
                    flow_rate_dispense = 100,
                    flow_rate_aspirate_mix = 25,
                    flow_rate_dispense_mix = 100,
                    air_gap_vol_bottom = 1,
                    air_gap_vol_top = 4,
                    disposal_volume = 1,
                    max_volume_allowed = 18,
                    reagent_volume = BEADS_VOLUME_PER_SAMPLE,
                    placed_in_multi = True,
                    first_well = 12,
                    v_fondo = 695) #1.95 * multi_well_rack_area / 2, #Prismatic
    
    Lysis = Simple_Reagent(name = 'Lysis',
                     flow_rate_aspirate     = 0.5,
                     flow_rate_dispense     = 0.5,
                     flow_rate_aspirate_mix = 0.5,
                     flow_rate_dispense_mix = 0.5,
                     delay                  = 0,
                     air_gap_vol_bottom     = 0
                     )

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('VALORES DE VARIABLES')
    ctx.comment(' ')
    ctx.comment('Número de muestras: ' + str(NUM_REAL_SAMPLES) + ' (' + str(num_cols) + ' columnas)')
    ctx.comment('Número de controles: ' + str(NUM_CONTROL_SPACES))
    ctx.comment(' ')
    ctx.comment('Número de mezclas en la muestra: ' + str(NUM_BEFORE_MIXES))
    ctx.comment('Número de mezclas en el deepwell: ' + str(NUM_AFTER_MIXES))
    ctx.comment(' ')
    ctx.comment('Volumen de muestra en el deepwell: ' + str(VOLUME_SAMPLE) + ' ul')
    ctx.comment('Volumen de solución con bolas magnéticas por muestra: ' + str(BEADS_VOLUME_PER_SAMPLE) + ' ul')
    ctx.comment('Volumen de proteinasa K por muestra: ' + str(PK_VOLUME_PER_SAMPLE) + ' ul')
    ctx.comment('Volumen de lysys por muestra: ' + str(LYSIS_VOLUME_PER_SAMPLE) + ' ul')
    ctx.comment(' ')
    ctx.comment('Número máximo de dispensaciones de lysys con la misma punta: ' + str(MAX_LYSYS_DISPENSE_PER_TIP) + ' muestras')
    ctx.comment(' ')
    ctx.comment('Foto-sensible: ' + str(PHOTOSENSITIVE))
    ctx.comment('Repeticiones del sonido final: ' + str(SOUND_NUM_PLAYS))
    ctx.comment(' ')

    
    ctx.comment('###############################################')
    ctx.comment('VOLUMENES PARA ' + str(num_samples) + ' muestras.')
    ctx.comment('')
    ctx.comment('Volumen de lysys necesario en B3 :' + str(LYSIS_VOLUME_PER_SAMPLE * NUM_REAL_SAMPLES) + " ul")
    Beads.comment_vol_info()
    Pk.comment_vol_info()
    ctx.comment('')

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

    def shake_pipet (pipet, rounds = 2, speed = 100, v_offset = 0):
        ctx.comment("Shaking " + str(rounds) + " rounds.")
        for i in range(rounds):
                pipet.touch_tip(speed = speed, radius = 0.1, v_offset = v_offset)

    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, drop_height, blow_out, touch_tip, x_offset,
                       pickup_height = 0, rinse = False, rinse_rounds = 2, mix_height = 0,
                       skipFinalAirGap = False, touch_tip_offset = -10, shake = False):
        '''
        x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        pickup_height: height from bottom where volume
        rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        drop_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        blow_out, touch_tip: if True they will be done after dispensing
        '''
        # Rinse before aspirating
        if rinse == True:
            custom_mix(pipet, reagent, location = source, vol = vol,
                       rounds = rinse_rounds, blow_out = True, mix_height = mix_height,
                       x_offset = x_offset)

        # SOURCE
        s = source.bottom(pickup_height).move(Point(x = x_offset[0]))
        if reagent.air_gap_vol_top > 0:
            ctx.comment ("Airgap top: " + str(reagent.air_gap_vol_top) + " ul")
            pipet.aspirate(reagent.air_gap_vol_top, source.top(z = -2),
                           rate = reagent.flow_rate_aspirate)  # air gap
        pipet.aspirate(vol, s, rate = reagent.flow_rate_aspirate)  # aspirate liquid
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pipet.aspirate(air_gap_vol, source.top(z = -2),
                           rate = reagent.flow_rate_aspirate)  # air gap

        # GO TO DESTINATION
        drop = dest.top(z = drop_height).move(Point(x = x_offset[1]))
        pipet.dispense(vol + air_gap_vol + reagent.air_gap_vol_top, drop,
                       rate = reagent.flow_rate_dispense)  # dispense all

        ctx.delay(seconds =reagent.delay) # pause for x seconds depending on reagent

        if shake == True:
            shake_pipet(pipet, rounds = 2, speed = 100, v_offset = drop_height)

        if blow_out == True:
            ctx.comment("Blowing out.") 
            pipet.blow_out(dest.top(z = drop_height))
        
        if touch_tip == True:
            pipet.touch_tip(speed = 20, radius = 0.8, v_offset = touch_tip_offset)

        if air_gap_vol != 0 and not skipFinalAirGap:
            #pipet.move_to(dest.top(z = drop_height))
            pipet.air_gap(air_gap_vol, height = drop_height) #air gap
        
        #if skipFinalAirGap and air_gap_vol != 0:
            #pipet.air_gap(air_gap_vol, height = drop_height) #air gap
            #pipet.move_to(dest.top(z = drop_height))
            #ctx.comment("Prueba de AirGap")
            #pipet.aspirate(air_gap_vol, dest.top(z = drop_height),rate = reagent.flow_rate_dispense)

    
    def distribute_custom(pipette, reagent, volume, src, dest, waste_pool, pickup_height, extra_dispensal, dest_x_offset, drop_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        pipette.aspirate((len(dest) * volume) +extra_dispensal
                         , src.bottom(pickup_height), rate = reagent.flow_rate_aspirate)
        pipette.move_to(src.top(z=5))
        #pipette.aspirate(reagent.air_gap_vol_bottom, rate = reagent.flow_rate_aspirate)  # air gap
        ctx.delay(seconds = 2) # pause for x seconds 
        pipette.air_gap(reagent.air_gap_vol_bottom, height = 5) #air gap
        for d in dest:
            pipette.dispense(volume + reagent.air_gap_vol_bottom, d.top(z = drop_height), rate = reagent.flow_rate_dispense)
            #pipette.move_to(d.top(z=5))
            #pipette.aspirate(reagent.air_gap_vol_bottom, rate = reagent.flow_rate_dispense)  # air gap
            pipette.air_gap(reagent.air_gap_vol_bottom, height = 15) #air gap
        try:
            pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
        except:
            pipette.blow_out(waste_pool.top(pickup_height + 3))
        return (len(dest) * volume)

    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 5, air_gap_vol = 0, blow_out_from_top = -15):
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
                        z = source_height).move(Point(x = x_offset[0])), rate = reagent.flow_rate_aspirate_mix)

        for _ in range(rounds):
            pipet.aspirate(vol, location = location.bottom(
                z = source_height).move(Point(x = x_offset[0])), rate = reagent.flow_rate_aspirate_mix)
            pipet.dispense(vol, location = location.bottom(
                z = mix_height).move(Point(x = x_offset[1])), rate = reagent.flow_rate_dispense_mix)

        pipet.dispense(1, location = location.bottom(
            z = mix_height).move(Point(x = x_offset[1])), rate = reagent.flow_rate_dispense_mix)

        if blow_out == True:
            pipet.blow_out(location.top(z = blow_out_from_top))  # Blow out
        
        if air_gap_vol > 0:
            pipet.air_gap(air_gap_vol, height = 0) #air gap

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
    
    def pick_up_tip(pip, tips = None):
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
    
    
    def calc_height(reagent, cross_section_area, aspirate_volume, min_height = 0.4 ):
        nonlocal ctx
        ctx.comment('¿Volumen útil restante ' + str(reagent.vol_well - reagent.dead_vol) +
                    ' uL < volumen necesario ' + str(aspirate_volume - reagent.disposal_volume * 8) + ' uL?')
        if (reagent.vol_well - reagent.dead_vol + 1) < (aspirate_volume - reagent.disposal_volume * 8):
            ctx.comment('Se debe utilizar el siguiente canal')
            ctx.comment('Canal anterior: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            ctx.comment('Nuevo canal: ' + str(reagent.col))
            reagent.vol_well = reagent.vol_well_original
            ctx.comment('Nuevo volumen: ' + str(reagent.vol_well) + ' uL')
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
            reagent.vol_well = reagent.vol_well - (aspirate_volume - reagent.disposal_volume * 8)
            ctx.comment('Volumen restante: ' + str(reagent.vol_well) + ' uL')
            if height < min_height:
                height = min_height
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
            reagent.vol_well = reagent.vol_well - (aspirate_volume - (reagent.disposal_volume * 8))
            ctx.comment('La altura calculada es ' + str(round(height, 2)) + ' mm')
            if height < min_height:
                height = min_height
            ctx.comment('La altura utilizada es ' + str(round(height, 2)) + ' mm')
            col_change = False
        return height, col_change
    
    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up_tip(pip, position = None):
        nonlocal tip_track
        #if not ctx.is_simulating():
        if recycle_tip:
            pip.pick_up_tip(pip.tip_racks[0].wells()[0])
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
                ctx.pause('Reemplaza las cajas de puntas de ' + str(pip.max_volume) + 'µl antes de continuar.')
                pip.reset_tipracks()
                tip_track['counts'][pip] = 0
                tip_track['num_refills'][pip] += 1
            if position is None:
                pip.pick_up_tip()
            else:
                pip.pick_up_tip(position)

    def drop_tip(pip, recycle = False, increment_count = True):
        nonlocal tip_track
        #if not ctx.is_simulating():
        if recycle or recycle_tip:
            pip.return_tip()
        else:
            pip.drop_tip(home_after = False)
        if increment_count:
            tip_track['counts'][pip] += 8


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
        '1000µl filter tiprack') for slot in ['10']]

    tips20 = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot, ' filter tiprack')
        for slot in ['11']]

    ################################################################################
    # setup samples and destinations
    sample_sources_full = generate_source_table(source_racks)
    sample_sources      = sample_sources_full[NUM_CONTROL_SPACES:num_samples]
    destinations        = dest_plate.wells()[NUM_CONTROL_SPACES:num_samples]
    destinations_full   = dest_plate.rows()[0][:num_samples]
    lysys_source        = lysys_rack.wells_by_name()['B3']
    dests_lysis         = list(divide_destinations(destinations, size_transfer))
    
    ctx.comment("***************************************************************************")
    ctx.comment("pk.fristwell: " + str(Pk.first_well))

    beads_reservoir = reagent_res.rows()[0][Beads.first_well - 1:Beads.first_well -1 + Beads.num_wells]
    pk_reservoir = reagent_res.rows()[0][Pk.first_well - 1:Pk.first_well - 1 + Pk.num_wells]

    p1000 = ctx.load_instrument(
        'p1000_single_gen2', 'right', 
        tip_racks = tips1000) # load P1000 pipette

    m20 = ctx.load_instrument(
        'p20_multi_gen2', 'left', 
        tip_racks = tips20) # load m20 pipette

    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {p1000: 0, m20: 0},
        'maxes': {p1000: 96 * len(p1000.tip_racks), m20: 96 * len(m20.tip_racks)}, #96 tips per tiprack * number or tipracks in the layout
        'num_refills' : {p1000 : 0, m20: 0},
        'tips': { p1000: [tip for rack in tips1000 for tip in rack.rows()[0]],
                    m20: [tip for rack in tips20 for tip in rack.rows()[0]]
                }

    }


    start_run()
    ############################################################################
    # STEP 1: ADD LYSIS 
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = log_step_start()

        used_vol = []
        dest_count = 0 # number of pipet uses

        for dest in dests_lysis:
            if not p1000.hw_pipette['has_tip']:
                 pick_up_tip(p1000)
            dest_count = dest_count + len(dest)

            num_mixes = 1
            #ctx.comment("Mezclas-   " + str(num_mixes))
            #custom_mix(p1000, reagent = Lysis, location = lysys_source, vol = LYSIS_VOLUME_PER_SAMPLE, 
            #    rounds = num_mixes, blow_out = False, mix_height = 15, x_offset = x_offset)

            used_vol_temp = distribute_custom(p1000, Lysis, volume = LYSIS_VOLUME_PER_SAMPLE,
                src = lysys_source, dest = dest,
                waste_pool = lysys_source, pickup_height = 2,
                extra_dispensal = extra_dispensal, dest_x_offset = 2, drop_height = 10)
            used_vol.append(used_vol_temp)

            # Check if it's time to change tip
            if dest_count >= MAX_LYSYS_DISPENSE_PER_TIP:
                ctx.comment("Changing tip, used " + str(dest_count) +" times.")
                drop_tip (p1000,recycle= recycle_tip)
                dest_count = 0

        drop_tip (p1000,recycle= recycle_tip)
        tip_track['counts'][p1000] += 1

        log_step_end(start)

    ############################################################################
    # STEP 2: MIX AND MOVE SAMPLES
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = log_step_start()

        for s, d in zip(sample_sources, destinations):
            if not p1000.hw_pipette['has_tip']:
                pick_up(p1000)

            # Mix the sample BEFORE dispensing
            if NUM_BEFORE_MIXES > 0:
                ctx.comment("Mezclas en origen " + str(NUM_BEFORE_MIXES))
                custom_mix(p1000, reagent = Samples, location = s, vol = volume_mix_tuberack, 
                    rounds = NUM_BEFORE_MIXES, blow_out = False, mix_height = 15, x_offset = x_offset)

            move_vol_multichannel(p1000, reagent = Samples, source = s, dest = d,
                vol = VOLUME_SAMPLE, air_gap_vol = Samples.air_gap_vol_bottom, x_offset = x_offset,
                pickup_height = 3, rinse = False, drop_height = -10,
                blow_out = NUM_AFTER_MIXES < 1, touch_tip = False, skipFinalAirGap = True)

            # Mix the sample BEFORE dispensing
            if NUM_AFTER_MIXES > 0:
                ctx.comment("Mezclas en destino " + str(NUM_AFTER_MIXES))
                custom_mix(p1000, reagent = Samples, location = d, vol = volume_mix_deepwell, 
                    rounds = NUM_AFTER_MIXES, blow_out = True, mix_height = 1,source_height = 1, x_offset = x_offset, air_gap_vol = 2 )

            p1000.drop_tip(home_after = False)
            tip_track['counts'][p1000] += 1

        log_step_end(start)

    
    ###############################################################################
    # STEP 3 Transferir Proteinasa K
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = log_step_start()

        pk_trips = math.ceil(Pk.reagent_volume / Pk.max_volume_allowed)
        pk_volume = Pk.reagent_volume / pk_trips
        pk_transfer_vol = []
        for i in range(pk_trips):
            pk_transfer_vol.append(pk_volume)
        
        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m20.hw_pipette['has_tip']:
                pick_up_tip(m20)
            for j,transfer_vol in enumerate(pk_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                # transfer_vol_extra = transfer_vol if j > 0 else transfer_vol + 100  # Extra 100 isopropanol for calcs
                # [pickup_height, change_col] = calc_height(Pk, multi_well_rack_area, transfer_vol_extra * 8)
                [pickup_height, change_col] = calc_height(Pk, multi_well_rack_area, transfer_vol * 8)
                
                ctx.comment('Aspirando desde la columna del reservorio: ' + str(Pk.first_well + Pk.col))
                ctx.comment('La altura de recogida es ' + str(round(pickup_height, 2)) + ' mm')
                move_vol_multichannel(m20, reagent = Pk, source = pk_reservoir[Pk.col],
                        dest = destinations_full[i], vol = transfer_vol, shake = True,
                        pickup_height = pickup_height, blow_out = True, touch_tip = False, drop_height = 5, 
                        air_gap_vol = Pk.air_gap_vol_bottom, x_offset = x_offset, skipFinalAirGap = True)

            # m20.air_gap(Pk.air_gap_vol_bottom, height = 5) #air gap

        if recycle_tip:
            m20.return_tip()
        else :
            drop_tip(m20)

        log_step_end(start)

    
    ###############################################################################
    # STEP 4 Transferir bolas magnéticas
    ###############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = log_step_start()

        beads_trips = math.ceil(Beads.reagent_volume / Beads.max_volume_allowed)
        beads_volume = Beads.reagent_volume / beads_trips
        beads_transfer_vol = []
        for i in range(beads_trips):
            beads_transfer_vol.append(beads_volume)
        
        rinse = True # Rinse first time
        rinse_rounds = 5
        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m20.hw_pipette['has_tip']:
                pick_up_tip(m20)

            for j,transfer_vol in enumerate(beads_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                # transfer_vol_extra = transfer_vol if j > 0 else transfer_vol + 100  # Extra 100 isopropanol for calcs
                # [pickup_height, change_col] = calc_height(Beads, multi_well_rack_area, transfer_vol_extra * 8)
                [pickup_height, change_col] = calc_height(Beads, multi_well_rack_area, transfer_vol * 8)
                
                ctx.comment('Aspirando desde la columna del reservorio: ' + str(Beads.first_well + Beads.col))
                ctx.comment('La altura de recogida es ' + str(round(pickup_height, 2)) + ' mm')
                move_vol_multichannel(m20, reagent = Beads, source = beads_reservoir[Beads.col],
                        dest = destinations_full[i], vol = transfer_vol, 
                        rinse = rinse, rinse_rounds = rinse_rounds, mix_height = 0,
                        touch_tip = False, touch_tip_offset = -20, shake = True,
                        pickup_height = pickup_height, blow_out = True, drop_height = 5, 
                        air_gap_vol = Beads.air_gap_vol_bottom, x_offset = x_offset, skipFinalAirGap = True)
                rinse_rounds = 1 # Disable Rinse

            # m20.air_gap(Beads.air_gap_vol_bottom, height = 5) #air gap

        if recycle_tip:
            m20.return_tip()
        else :
            drop_tip(m20)

        log_step_end(start)

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

    finish_run(switch_off_lights)