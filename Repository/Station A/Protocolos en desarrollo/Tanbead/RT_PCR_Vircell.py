from opentrons import protocol_api

metadata = {
    'protocolName': 'Multidispensacion',
    'author': 'RGC',
    'description': 'Preparación placa Tan Bead',
    'apiLevel': '2.7'
    
    }
samples1=1
samples2=4
samples3=7
samples4=10
samples5=11
k=6
samples6=8
tips1=9
placa=3
samples7=5
tips2=2


def run(protocol: protocol_api.ProtocolContext):
#labware
    deepwell = protocol.load_labware('vircell_96_wellplate_2000ul', placa)
    proteink = protocol.load_labware('vircell_15_tuberack_2000ul', k)
    samplesplate1 = protocol.load_labware('vircell_15_tuberack_15000ul', samples1)
    samplesplate2 = protocol.load_labware('vircell_15_tuberack_15000ul', samples2)
    samplesplate3 = protocol.load_labware('vircell_15_tuberack_15000ul', samples3)
    samplesplate4 = protocol.load_labware('vircell_15_tuberack_15000ul', samples4)
    samplesplate5 = protocol.load_labware('vircell_15_tuberack_15000ul', samples5)
    samplesplate6 = protocol.load_labware('vircell_15_tuberack_15000ul', samples6)
    samplesplate7= protocol.load_labware('vircell_15_tuberack_15000ul', samples7)
    tiprack_1 = protocol.load_labware('Opentrons_96_tiprack_300ul', tips1)
    tiprack_2 = protocol.load_labware('opentrons_96_filtertiprack_1000ul', tips2)

#pipettes
    right_pipette = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack_1])
    left_pipette = protocol.load_instrument('p1000_single_gen2', 'left', tip_racks=[tiprack_2])

#commands

    protocol.set_rail_lights(True)
    
    right_pipette.distribute(
            10,
            [proteink.wells_by_name()['A1']],[deepwell.columns_by_name() [x] for x in ['1','2','3','4','5','6','7','8','9','10','11','12']],
            disposal_volume=3
            )
            
    left_pipette.transfer(
            300,
            [samplesplate1.wells_by_name() [x] for x in ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5']],
            [deepwell.wells_by_name() [y] for y in ['A1','B1','C1','D1','E1','F1','G1','H1','A2','B2','C2','D2','E2','F2','G2']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )
    left_pipette.transfer(
            300,
            [samplesplate2.wells_by_name() [x] for x in ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5']],
            [deepwell.wells_by_name() [y] for y in ['H2','A3','B3','C3','D3','E3','F3','G3','H3','A4','B4','C4','D4','E4','F4']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )
    left_pipette.transfer(
            300,
            [samplesplate3.wells_by_name() [x] for x in ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5']],
            [deepwell.wells_by_name() [y] for y in ['G4','H4','A5','B5','C5','D5','E5','F5','G5','H5','A6','B6','C6','D6','E6']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )
    left_pipette.transfer(
            300,
            [samplesplate4.wells_by_name() [x] for x in ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5']],
            [deepwell.wells_by_name() [y] for y in ['F6','G6','H6','A7','B7','C7','D7','E7','F7','G7','H7','A8','B8','C8','D8']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )
    left_pipette.transfer(
            300,
            [samplesplate5.wells_by_name() [x] for x in ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5']],
            [deepwell.wells_by_name() [y] for y in ['E8','F8','G8','H8','A9','B9','C9','D9','E9','F9','G9','H9','A10','B10','C10']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )
    left_pipette.transfer(
            300,
            [samplesplate6.wells_by_name() [x] for x in ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5']],
            [deepwell.wells_by_name() [y] for y in ['D10','E10','F10','G10','H10','A11','B11','C11','D11','E11','F11','G11','H11','A12','B12']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )
    left_pipette.transfer(
            300,
            [samplesplate7.wells_by_name() [x] for x in ['A1','A2','A3','A4']],
            [deepwell.wells_by_name() [y] for y in ['C12','D12','E12','F12']],
            disposal_volume=0,
            blow_out=True,
            new_tip='always',
            air_gap=2
            )




















     
