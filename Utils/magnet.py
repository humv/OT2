from opentrons import protocol_api

metadata = {'apiLevel': '2.3'}

def run(protocol: protocol_api.ProtocolContext):
    mag_mod = protocol.load_module('Magnetic Module Gen2', '1')
    plate = mag_mod.load_labware('kingfisher_96_wellplate_2000ul', 'KingFisher 96 Well Plate 2mL')

    for i in range(50):
        protocol.pause()
        mag_mod.engage(height = 7)
        protocol.pause()
        mag_mod.disengage()