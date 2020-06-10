from opentrons import protocol_api

# Metadatos:
metadata = {
    'protocolName': 'Protocolo Extracción Turbobeads',
    'apiLevel': '2.2'}


def run(protocol: protocol_api.ProtocolContext):
    tiprack1 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '1')
    tiprack2 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '2')
    tiprack3 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '3')
    tiprack4 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '4')
    tiprack5 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '5')
    tiprack6 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '6')
    reservoir7 = protocol.load_labware('nest_12_reservoir_15ml', '7')
    reservoir8 = protocol.load_labware('nest_1_reservoir_195ml', '8')
    wellplate9 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '9')
    mag_mod = protocol.load_module('magnetic module', '10')
    wellplate10 = mag_mod.load_labware('usascientific_96_wellplate_2.4ml_deep')
    reservoir11 = protocol.load_labware('nest_1_reservoir_195ml', '11')

    p300multi = protocol.load_instrument('p300_multi_gen2', 'right', tip_racks=[tiprack4, tiprack5, tiprack6, tiprack3, tiprack2, tiprack1])

    comentario = '----------> 0.'
    protocol.comment(comentario)
    mag_mod.disengage()
    p300multi.well_bottom_clearance.aspirate = 1
    p300multi.well_bottom_clearance.dispense = 30
    blb_lane = reservoir7.columns_by_name()['2']

    p300multi.pick_up_tip()
    for blb in range(12):
        if blb ==4:
            blb_lane = reservoir7.columns_by_name()['3']
        elif blb == 8:
            blb_lane = reservoir7.columns_by_name()['4']
        for _ in range(2):
            p300multi.transfer(150, blb_lane, wellplate10.columns()[blb], new_tip='never')
            p300multi.blow_out(wellplate10.wells()[blb * 8].top(z=-15))
    p300multi.drop_tip()


    comentario = '----------> 1.'
    protocol.comment(comentario)
    p300multi.pick_up_tip()
    for bead in range(12):
        p300multi.transfer(20, reservoir7.columns_by_name()['10'], wellplate10.columns()[bead], new_tip='never')
        p300multi.blow_out(wellplate10.wells()[bead * 8].top(z=-15))
    p300multi.drop_tip()


    comentario = '----------> 2.'
    protocol.comment(comentario)
    isop_lane = reservoir7.columns_by_name()['6']

    for isop in range(12):
        p300multi.pick_up_tip()
        if isop ==4:
            isop_lane = reservoir7.columns_by_name()['7']
        elif isop == 8:
            isop_lane = reservoir7.columns_by_name()['8']
        p300multi.well_bottom_clearance.dispense = 30
        p300multi.transfer(200, isop_lane, wellplate10.columns()[isop], new_tip='never')
        p300multi.blow_out(wellplate10.wells()[isop * 8].top(z=-15))

        p300multi.well_bottom_clearance.dispense = 5
        p300multi.transfer(200, isop_lane, wellplate10.columns()[isop], new_tip='never', mix_after=(12, 180))
        p300multi.blow_out(wellplate10.wells()[isop * 8].top(z=-15))
        p300multi.drop_tip()


    comentario = '----------> 3.'
    protocol.comment(comentario)
    protocol.delay(minutes=5)


    comentario = '----------> 4.'
    protocol.comment(comentario)
    mag_mod.engage() ############ ALTURA
    protocol.delay(minutes=3)
    p300multi.well_bottom_clearance.aspirate = 0
    p300multi.well_bottom_clearance.dispense = 20

    for sobre1 in range(12):
        p300multi.pick_up_tip()
        for _ in range(5):
            p300multi.transfer(200, wellplate10.columns()[sobre1], reservoir11.columns_by_name()['1'], new_tip='never')
            p300multi.blow_out(reservoir11.wells()[0].top(z=-15))
        p300multi.drop_tip()


    comentario = '----------> 5. '
    protocol.comment(comentario)
    mag_mod.disengage()
    p300multi.well_bottom_clearance.aspirate = 1

    for wf1 in range(12):
        p300multi.pick_up_tip()
        p300multi.well_bottom_clearance.dispense = 30
        p300multi.transfer(150, reservoir8.columns_by_name()['1'], wellplate10.columns()[wf1], new_tip='never')
        p300multi.blow_out(wellplate10.wells()[wf1 * 8].top(z=-15))

        p300multi.well_bottom_clearance.dispense = 5
        p300multi.transfer(150, reservoir8.columns_by_name()['1'], wellplate10.columns()[wf1], new_tip='never', mix_after=(10, 180))
        p300multi.blow_out(wellplate10.wells()[wf1 * 8].top(z=-15))
        p300multi.drop_tip()


    comentario = '----------> 6.'
    protocol.comment(comentario)
    mag_mod.engage() ############ ALTURA
    p300multi.well_bottom_clearance.aspirate = 0
    p300multi.well_bottom_clearance.dispense = 20

    for sobre2 in range(12):
        p300multi.pick_up_tip()
        for _ in range(2):
            p300multi.transfer(160, wellplate10.columns()[sobre2], reservoir11.columns_by_name()['1'], new_tip='never')
            p300multi.blow_out(reservoir11.wells()[0].top(z=-15))
        p300multi.drop_tip()


    comentario = '----------> 7.'
    protocol.comment(comentario)
    mag_mod.disengage()
    p300multi.well_bottom_clearance.aspirate = 1

    for wf2 in range(12):
        p300multi.pick_up_tip()
        p300multi.well_bottom_clearance.dispense = 30
        p300multi.transfer(150, reservoir8.columns_by_name()['1'], wellplate10.columns()[wf2], new_tip='never')
        p300multi.blow_out(wellplate10.wells()[wf2 * 8].top(z=-15))

        p300multi.well_bottom_clearance.dispense = 5
        p300multi.transfer(150, reservoir8.columns_by_name()['1'], wellplate10.columns()[wf2], new_tip='never', mix_after=(10, 180))
        p300multi.blow_out(wellplate10.wells()[wf2 * 8].top(z=-15))
        p300multi.drop_tip()




    protocol.pause('VOLVER A PONER CAJAS DE PUNTAS LLENAS EN SLOTS 4, 5 y 6')
    p300multi.reset_tipracks()




    comentario = '----------> 8.'
    protocol.comment(comentario)
    mag_mod.engage() ############ ALTURA
    p300multi.well_bottom_clearance.aspirate = 0
    p300multi.well_bottom_clearance.dispense = 20

    for sobre2 in range(12):
        p300multi.pick_up_tip()
        for _ in range(2):
            p300multi.transfer(160, wellplate10.columns()[sobre2], reservoir11.columns_by_name()['1'], new_tip='never')
            p300multi.blow_out(reservoir11.wells()[0].top(z=-15))
        p300multi.drop_tip()


    comentario = '----------> 9.'
    protocol.comment(comentario)
    protocol.delay(minutes=2)


    comentario = '----------> 10.'
    protocol.comment(comentario)
    p300multi.well_bottom_clearance.aspirate = 1
    p300multi.well_bottom_clearance.dispense = 1

    for elution in range(12):
        mag_mod.disengage()
        p300multi.pick_up_tip()
        p300multi.transfer(40, reservoir7.columns_by_name()['12'], wellplate10.columns()[elution], new_tip='never', mix_after=(10, 30))
        p300multi.blow_out(wellplate10.wells()[elution * 8].top(z=-15))
        p300multi.drop_tip()


    comentario = '----------> 11.'
    protocol.comment(comentario)
    mag_mod.engage()  ############ ALTURA
    protocol.delay(minutes=2)

    for elution2 in range(12):
        p300multi.pick_up_tip()
        p300multi.well_bottom_clearance.dispense = 5
        p300multi.transfer(40, wellplate10.columns()[elution2], wellplate9.columns()[elution2], new_tip='never')
        p300multi.blow_out(wellplate9.wells()[elution2 * 8].top(z=-5))
        p300multi.drop_tip()