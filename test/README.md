## Test directory

### Run unit tests and intergration test

Run the following command to test:

```sh
~$ make
```

It will run the unit tests and integration tests if they are available.
Run 'make help' for more details and options.

### Generate GDB version of the programs for debug

To generate a GDB version of the programs.
Run the following make command:

```sg
~$ make progname
```

The program make will generate a GDB version of the program and load the arguments using a GDB script.

You are ready to debug sfstereoniptomo with GDB.
