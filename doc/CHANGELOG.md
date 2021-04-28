# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## UNRELEASED - Contribute with the next version

## [v0.1.3-alpha.1](https://github.com/Dirack/creVelocityInversion/compare/v0.1.2-alpha.1...develop/0.1.3) - 2021-04-28

## RELEASED

## [v0.1.2-alpha.1](https://github.com/Dirack/creVelocityInversion/releases/tag/v0.1.2-alpha.1) (Development) - 2021-04-28

## Added

- Velocity funtion to update the velocity gradient to simulate a velocity varying with depth model
and evaluate if the NIP sources converges to the reflector interfaces [2c306e3](https://github.com/Dirack/creVelocityInversion/commit/2c306e3)

## Changed

- Run inversion in looping [#20](https://github.com/Dirack/creVelocityInversion/issues/20) Run the velocity inversion
in a loop using the prevous optimal model as input for the next iteration [87a0567](https://github.com/Dirack/creVelocityInversion/commit/87a0567)

## Removed

- Remove program sfniptimecurve and sfgetbetaangle: Those were test programs and there is no need for them anymore
[8b0b723](https://github.com/Dirack/creVelocityInversion/commit/8b0b723)

- Remove useless experiments from experiments directory: The inversion strategy has changed from diffraction simulation
to ray tracying algorithm and NIP tomography [be73a35](https://github.com/Dirack/creVelocityInversion/commit/be73a35)

## [v0.1.1-alpha.1](https://github.com/Dirack/creVelocityInversion/releases/tag/v0.1.1-alpha.1) (Development) - 2020-10-07

## Added

- Scratch of the velocity inversion algorithm with iterative picking [#2](https://github.com/Dirack/creVelocityInversion/issues/2) [42c68d5](https://github.com/Dirack/creVelocityInversion/commit/42c68d5)

- Use the Cameron (2008) velocity inversion algorithm Resolve [#9](https://github.com/Dirack/creVelocityInversion/issues/9) [3d7331e](https://github.com/Dirack/creVelocityInversion/commit/3d7331e) [2800abb](https://github.com/Dirack/creVelocityInversion/commit/2800abb)

- Add parameters dictionary to allow test and command line input parameters [0e28fc0](https://github.com/Dirack/creVelocityInversion/commit/0e28fc0)

- Casting parameters function [5090942](https://github.com/Dirack/creVelocityInversion/commit/5090942)
