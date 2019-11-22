## Working with data

### Transferring data

You can use GUI software to transfer data files:
- For transferring data from your computer to a remote storage location, we recommend [Cyberduck](https://sciwiki.fredhutch.org/compdemos/Mountain-CyberDuck/#cyberduck)
- For transferring data among Fred Hutch storage locations, we recommend [Motuz](https://sciwiki.fredhutch.org/scicomputing/store_overview/#data-transfer)

Of course, there are also Unix commands to transfer files, called `scp` (secure copy).

  scp

You can copy entire directories by making your command recursive:

  scp -r

### Working with large and/or many data files

* uncompressing zip files from the command line

  zip
  gunzip

  tar

## Environments

### Environment variables

What is an environment?

* environment variables; PATH

### Finding and loading software

  module load
  ml


  module spider

## Running software

* processes: foreground, background, killing, `top`
* grabnode and cluster job submission


* ssh access keys
