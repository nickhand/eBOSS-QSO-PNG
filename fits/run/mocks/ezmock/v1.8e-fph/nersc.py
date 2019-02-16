from argparse import ArgumentParser
from nbodykit.lab import TaskManager
from eboss_qso.fits.driver import QSOFitDriver
import importlib, os

def main(ns):

    # the module name
    modname = os.path.splitext(os.path.split(ns.filename)[-1])[0]
    mod = importlib.import_module(modname)
    add_commands = getattr(mod, 'add_commands')
    RSDFitRunner = getattr(mod, 'RSDFitRunner')

    # run with 1 core per task
    with TaskManager(cpus_per_task=1, use_all_cpus=True) as tm:

        # whether we are varying shot noise
        vary_shot_noise = False if ns.fix_shot_noise else True

        if len(ns.files):
            files = ns.files
        else:
            files = range(ns.start, ns.stop, ns.step)

        # iterate through the boxes in parallel
        for box in tm.iterate(files):

            # clear old commands and update the list of commands for this box
            RSDFitRunner.commands.clear()
            add_commands(box=box, vary_shot_noise=vary_shot_noise, kmin=ns.kmin,
                            cov=ns.cov, kmax=ns.kmax, use_temp_files=True)

            # get the command and run
            command = RSDFitRunner.commands[ns.testno]
            print("executing %s..." %command)
            QSOFitDriver.run_from_args(command.split()[1:], comm=tm.comm)


if __name__ == '__main__':

    from eboss_qso.fits import NERSCManager

    h = 'the name of the file to run'
    NERSCManager.add_argument('filename', type=str, help=h)
    NERSCManager.add_argument('testno', type=int, help=h)

    NERSCManager.add_argument('--start', type=int, required=True)
    NERSCManager.add_argument('--stop', type=int, required=True)
    NERSCManager.add_argument('--step', type=int, default=1)

    NERSCManager.add_argument('--files', nargs="*", type=int, default=[])

    # other config
    NERSCManager.add_argument('--cov', choices=['mock', 'analytic'], required=True)
    NERSCManager.add_argument('--fix-shot-noise', action='store_true', default=False)
    NERSCManager.add_argument('--kmin', type=float, default=1e-4)
    NERSCManager.add_argument('--kmax', type=float, default=0.3)

    ns = NERSCManager.parse_args()
    main(ns)
