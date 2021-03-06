This is a experiment for sfstereoniptomo test with GDB.
Run this experiment before test to generate the input.

After that, go to the TDD directory, build the program version for
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
