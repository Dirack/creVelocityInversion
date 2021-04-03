## Test directory

### Run unit tests and intergration test

Run the following command to test:

```sh
~$ make
```

It will run the unit tests of diffmig and diffsimul experiments using unittest python module and it will run optimization tests to get the time of those experiments.

### Generate GDB version of the programs for debug

To generate a GDB version of the programs.
Run the following make command:

```sg
~$ make progname.x
```

The program make will generate a GDB version of the program and load the arguments using a GDB script.

You are ready to debug sfstereoniptomo with GDB.
