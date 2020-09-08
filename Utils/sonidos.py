import opentrons.execute
import subprocess
from opentrons import protocol_api

AUDIO_FILE_PATH1 = '/var/lib/jupyter/notebooks/sonidos/finished_process_esp.mp3'
AUDIO_FILE_PATH2 = '/var/lib/jupyter/notebooks/sonidos/finalizado.mp3'
def run_quiet_process(command):
     subprocess.check_output('{} &> /dev/null'.format(command), shell=True)
def test_speaker():
     print('Speaker')
     print('Next\t--> CTRL-C')
     try:
         run_quiet_process('mpg123 {}'.format(AUDIO_FILE_PATH1))
         run_quiet_process('mpg123 {}'.format(AUDIO_FILE_PATH2))
         run_quiet_process('mpg123 {}'.format(AUDIO_FILE_PATH1))
     except KeyboardInterrupt:
         pass
         print()

protocol = opentrons.execute.get_protocol_api('2.6')
protocol.home() 
test_speaker()
