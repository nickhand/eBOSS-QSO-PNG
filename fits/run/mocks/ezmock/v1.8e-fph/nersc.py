from argparse import ArgumentParser
from nbodykit.lab import TaskManager
from eboss_qso.fits.driver import QSOFitDriver
from run_fnl_fits import RSDFitRunner, add_commands

def main(ns):

    with TaskManager(cpus_per_task=1, use_all_cpus=True) as tm:

        vary_shot_noise = False if ns.fix_shot_noise else True

        for box in tm.iterate(range(ns.start, ns.stop, ns.step)):

            RSDFitRunner.commands.clear()
            add_commands(box=box, vary_shot_noise=vary_shot_noise, kmax=ns.kmax, use_temp_files=True)

            command = RSDFitRunner.commands[ns.testno]
            QSOFitDriver.run_from_args(command.split()[1:], comm=tm.comm)


if __name__ == '__main__':

    from eboss_qso.fits import NERSCManager

    h = 'the name of the file to run'
    NERSCManager.add_argument('filename', type=str, help=h)
    NERSCManager.add_argument('testno', type=int, help=h)

    NERSCManager.add_argument('--start', type=int, required=True)
    NERSCManager.add_argument('--stop', type=int, required=True)
    NERSCManager.add_argument('--step', type=int, default=1)
    NERSCManager.add_argument('--fix-shot-noise', action='store_true', default=False)
    NERSCManager.add_argument('--kmax', type=float, default=0.3)

    ns = NERSCManager.parse_args()
    main(ns)
