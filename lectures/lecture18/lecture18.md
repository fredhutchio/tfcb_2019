## Working with remote data files

### Transferring data

So far, we've used GitHub as a location to store files so we can access them remotely (via `git pull`). This has worked for the small, public files we've been using for class, but this is insufficient for most research projects.

You can use GUI software to transfer data files:
- For transferring data from your computer to a remote storage location, we recommend [Cyberduck](https://sciwiki.fredhutch.org/compdemos/Mountain-CyberDuck/#cyberduck)
- For transferring data among Fred Hutch storage locations, we recommend [Motuz](https://sciwiki.fredhutch.org/scicomputing/store_overview/#data-transfer)

Of course, there are also Unix commands to transfer files, called `scp` (secure copy). More complete instructions can be found [here](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/), but the most straightforward for our purposes will be to open our shell and transfer a file from our local computer (laptop) to the remote computer (rhino). Navigate to  `tcfb_2019/lectures/lecture18` on your local computer and then execute:

  scp setup.sh USERNAME@rhino:~

This places the file `setup.sh` in your home directory (`~`) on rhino.

You can copy entire directories by making your command recursive:

  scp -r PATH/TO/DIRECTORY USERNAME@rhino:PATH/TO/RHINO/DIRECTORY

  wget
  curl

### Working with large and/or many data files

* uncompressing zip files from the command line

  tar

  zip
  gunzip

## Environments

### Environment variables

What is an environment?

* environment variables; PATH

### Finding and loading software

  module list

  module spider

  module load
  ml

## Running software

* processes: foreground, background, killing, `top`
* grabnode and cluster job submission


* ssh access keys
