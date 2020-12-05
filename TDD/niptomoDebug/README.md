This is a experiment for sfniptomo test with GDB.
Run this experiment before test to generate the input.

After that, go to the TDD directory, build the program version for
GDB with:

```sg
~$ make niptomo.x
```

Call this program in GDB with:

```sh
~$ gdb niptomo.x
```

Use the gdb\_script inside gdb session to set up the arguments
of the program:

```
(gdb) source gdb_script
```

You are ready to debug sfniptomo with GDB.
