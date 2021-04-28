## Stereoniptomo test 

### Picking

To run this experiment, first do the time picking of the reflectors
in the stacked section using the program _sfipick_
with the following command:

```sh
~$ make picking
```

This command will convert the stacked section in ascii format to RSF
format and it will run the program _sfipick_ to do iterative picking. 

### Run full experiment

The file 'pick.txt' will be generated to store time picked points
from stacked section. After that, run full experiment with the command:

```sh
~$ make
```

That command will run the experiment using the points you picked in the
stacked section as a parameter. Those points are associated with the RNIP
and BETA parameters in the 'parameterscube.rsf' file that will be used to
set up the model tracing rays to locate the NIP sources using _sfnipmodsetup_.

The program _sfstereoniptomo_ does the global optimization of the NIP sources
location throushout forward modeling by a ray-tracing from NIP sources
to acquisition surface and using traveltime as convergence criteria.

### Generate a GDB version of sfstereoniptomo

To generate a sfstereoniptomo for GDB,
go to the TDD directory, build the program version for
GDB with:

```sg
~$ make stereoniptomo.x
```

Call this program in GDB with:

```sh
~$ gdb stereoniptomo.x
```

Use the gdb\_script inside gdb session to set up the arguments
of the program:

```
(gdb) source gdb_script5
```

You are ready to debug sfstereoniptomo with GDB.
