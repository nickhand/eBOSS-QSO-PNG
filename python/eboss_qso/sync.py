import argparse
import os

RSYNC = "rsync -e ssh -avzl --progress --delete"

class NERSCSync(object):
    """
    Sync to NERSC

    Parameters
    ----------
    direction : {'to', 'from'}
        the direction to do the transfer, either 'to' or 'from'
        the remote host
    host : {'cori', 'edison'}
        the name of the remote host
    path : str
        a subpath from the ``RSDFIT_FITS`` directory
    dry_run : bool, optional
        whether to do a dry-run
    delete : bool, optional
        whether to delete files when syncing
    """

    def __init__(self, direction, host, path=None, dry_run=False, delete=False):
        self.direction = direction
        self.host = host
        self.path = path
        self.dry_run = dry_run
        self.delete = delete

    def __call__(self, remote_dir, local_dir, exclude=[]):
        """
        Parameters
        ----------
        remote_dir : str
            the directory on NERSC to sync
        local_dir : str
            the local directory to to sync
        exclude : list of str, optional
            additional things to exclude
        """
        # the command + options
        if not self.delete:
            cmd = [RSYNC.split('--delete')[0]]
        else:
            cmd = [RSYNC]

        cmd += exclude
        if self.dry_run: cmd.append('--dry-run')

        # add the directories and run the command
        dirs = [local_dir, remote_dir]
        if self.path is not None:
            dirs[0] += self.path
            dirs[1] += self.path

            # append a backslash
            for i, d in enumerate(dirs):
                dirs[i] = os.path.abspath(d) + os.sep

        dirs[1] = "nhand@%s:%s" %(self.host, dirs[1])
        if self.direction == 'from':
            dirs = dirs[::-1]
        cmd += dirs
        ret = os.system(" ".join(cmd))

def sync_eboss():

    # the main parser
    desc = "sync the eBOSS results `rsync`"
    parser = argparse.ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(dest='subparser_name')

    # setup the parent
    parent = argparse.ArgumentParser(add_help=False)
    h = 'the transfer direction; either `to` or `from` the remote host'
    parent.add_argument('direction', type=str, choices=['to', 'from'], help=h)
    h = 'the remote host; either `cori` or `edison`'
    parent.add_argument('host', type=str, choices=['cori', 'edison'], help=h)
    h = 'an additional subpath to sync only'
    parent.add_argument('-d', '--dir', type=str, help=h)
    h = 'show what would have been transferred'
    parent.add_argument('-n', '--dry-run', action='store_true', help=h)
    h = 'add the --delete tag'
    parent.add_argument('--delete', action='store_true', help=h)

    # measurements
    h = 'sync the `measurements` directory'
    a = subparsers.add_parser('measurements', parents=[parent], help=h)

    # measurements
    h = 'sync the `fits` directory'
    a = subparsers.add_parser('fits', parents=[parent], help=h)

    # reports
    h = 'sync the `reports` directory'
    b = subparsers.add_parser('reports', parents=[parent], help=h)

    # notebooks
    h = 'sync the `notebooks` directory'
    b = subparsers.add_parser('notebooks', parents=[parent], help=h)

    ns = parser.parse_args()
    sync = NERSCSync(ns.direction, ns.host, path=ns.dir, dry_run=ns.dry_run, delete=ns.delete)

    if ns.subparser_name == 'measurements':
        remote_dir = "/global/cscratch1/sd/nhand/Research/eBOSS/measurements/"
        local_dir = "/Users/nhand/Research/Analysis/thesis/eBOSS-QSO-PNG/measurements/"
        exclude = ["--exclude='info'", "--exclude='plots'", "--exclude='.*'", '--exclude=run']
    elif ns.subparser_name == 'reports':
        remote_dir = "/global/project/projectdirs/m779/www/nhand/notebooks/eboss-qso-fits/"
        local_dir = "/Users/nhand/Research/Analysis/thesis/eBOSS-QSO-PNG/fits/reports/"
        exclude = ["--exclude='.*'"]
    elif ns.subparser_name == 'notebooks':
        remote_dir = "/global/project/projectdirs/m779/www/nhand/notebooks/eboss-qso-notebooks/"
        local_dir = "/Users/nhand/Research/Analysis/thesis/eBOSS-QSO-PNG/notebooks/"
        exclude = ["--exclude='.*'"]
    elif ns.subparser_name == 'fits':
        remote_dir = "/global/project/projectdirs/m779/nhand/Research/eBOSS/fits/results/"
        local_dir = "/Users/nhand/Research/Analysis/thesis/eBOSS-QSO-PNG/fits/results/"
        exclude = ["--exclude='.*'"]

    sync(remote_dir, local_dir, exclude=exclude)
