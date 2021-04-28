# CRE Velocity Inversion

> Velocity inversion algorithm using Common Reflection Element (CRE) traveltime approximation.

[Designed to Madagascar package](https://ahay.org)

[![Github release](https://img.shields.io/github/v/release/Dirack/creVelocityInversion)](https://github.com/Dirack/creVelocityInversion/releases/latest) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Madagascar](https://img.shields.io/badge/Madagascar-v3.0-blue)](https://github.com/ahay/src/tree/master) [![Build Status](https://travis-ci.org/Dirack/creVelocityInversion.svg?branch=master)](https://travis-ci.org/Dirack/creVelocityInversion)

Common Reflection Element (CRE) velocity inversion algorithm is based on [Cameron (2008)](http://www.reproducibility.org/RSF/book/tccs/time2depth/paper_html/) velocity inversion algorithm and 
in diffraction focusing migration and simulation. The depth velocity model is obtained from stacked section in two steps:
First, simulation and migration of diffraction hyperbolas in the stacked section to obtain the time migration velocity model.
Second, time to depth conversion of the time velocity model.

## Development setup

- Madagascar package (3.0)

You need to have the actual Madagascar package stable release installed on your computer. Please follow the
[Installing Madagascar page](http://www.ahay.org/wiki/Installation) in the official documentation. You can install
Madagascar automatically from Shell Script using program _madagainstall_ from [Shellinclude library](https://github.com/Dirack/Shellinclude/tree/v1.2.2-beta.1).

## Installation

After Madagascar installing process, you need to install the programs of this repository in your local Madagascar user's
directory. You can compile and install it as any other Madagascar program. 
Usually, Madagascar keeps the path of your local copy source files in the $RSFSRC environment variable. You can
show that on a bash terminal using 'echo' command:

```sh
~$ echo "$RSFSRC"
```

And Madagascar will install executable files on your $RSFROOT directory. You can show that environment variable
with 'echo' too:

```sh
~$ echo "$RSFROOT"
```

Madagascar stores user programs in $RSFSRC/user directory. So, you can create a new directory or put this
repository inside that directory. In this repository, such as every user's repository in Madagascar, we have a compilation 
[SConstruct](https://github.com/Dirack/vfsa/blob/master/SConstruct) that compile the C programs.
Run 'scons' on your $RSFSRC/user/creGatherInterpolation repository to compile it:

```shell
~$ scons
```

And run 'scons install' in the top directory of your local Madagascar installation 
(the directory path in your $RSFSRC variable):

```shell
~$ sudo scons install
```

If you have any doubt about this process, please reffer to the Madagascar oficial documentation in 
[Adding_new_programs_to_Madagascar](http://www.ahay.org/wiki/Adding_new_programs_to_Madagascar)

## Usage example

A few motivating and useful examples of how that product can be used. 
_For more examples and details, please refer to the [Wiki](https://github.com/Dirack/creVelocityInversion/wiki)._

We also have many SConstruct examples in this repository in the
[experiments directory](https://github.com/Dirack/creVelocityInversion/tree/master/experiments)

## Release History and Versions

This package version is referenced in VERSION.md file, and you can see the [complete release history in our wiki](https://github.com/Dirack/creVelocityInversion/wiki/Release-history) or in the CHANGELOG.md file.

## Meta

[main page](https://github.com/Dirack/creVelocityInversion)

Rodolfo Dirack – [@dirack](https://github.com/Dirack) – rodolfo_profissional@hotmail.com

Distributed under the GPL3 license. See ``LICENSE`` for more information.

## Contributing

In order to contribute with this project you should follow the list of steps bellow, please check out ["How to contribute with this project?"](https://github.com/Dirack/creVelocityInversion/wiki/Contribute) in our Wiki for more details. 

1. Create an issue to your request or choose an issue already defined
2. Fork this project in https://github.com/Dirack/creVelocityInversion/fork 
3. Create a branch for your contribution (name it using gitflow)
4. Do clear _commit_ messages (a title with 50 characters and two paragraphs of text)
5. _Push_ your contribution to this repository
6. Create a new Pull Request with a clear description of your contribution

###### Important: The commit history should be clear, with commit mesages around one or two paraghraps describing your modifications. Pull Requests with unsatisfactory commit history will be rejected.
