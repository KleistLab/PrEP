# prophylaxis estimator
>compute the prophylactic efficacy of DTG using different methods

## Table of Contents
- [CLI](#cli)
  * [extrande](#extrande)
  * [ntm](#ntm)
  * [ctsm](#ctsm)
  * [pgs](#pgs)  
- [Highlights](#highlights)
- [About](#about)

## CLI
```
usage: ./prophylaxisEstimator.py <command> [options]

command:
  Choose one method for computing

  {extrande,ntm,ctsm,pgs}
    extrande            run method extrande
    ntm                 run next transition method
    ctsm                run constant time step method
    pgs                 run probability generating system

Run './prophylaxisEstimator.py COMMAND --help' for more information on a command
```
### extrande
```
usage: ./prophylaxisEstimator.py extrande [options]

optional arguments:
  -h, --help        show this help message and exit
  --mdose float     mass of DTG dose [mg] (default: 50)
  --adh float       adherence of drug taken (default: 1)
  --cdose int       count of doses in a regimen (default: 3)
  --inputfile str   specify the inputfile that contains the pharmacokinetic
                    parameters. For deterministic methods please input only                           
                    one set of pk parameters (default: PK_default.csv)                                
  --outputfile str  specify the name of output csv file (default:
                    extrande.csv)
  --ifpdf bool      if output the graphic of results in pdf (default: False)
  --csimul int      count of simulations (default: 5000)
  --tps list        time points to begin simulation, i.e. time of viral
                    exposure relative to timing of first dose [hr] (default:
                    [1,6,18,23])
  --vload str       int or str {"homo", "hetero"}, initial viral load, can be
                    fixed number or generated randomly according to the
                    transition mode (default: 1)
```
### ntm
```
usage: ./prophylaxisEstimator.py ntm [<options>]

optional arguments:
  -h, --help            show this help message and exit
  --mdose float         mass of DTG dose [mg] (default: 50)
  --adh float           adherence of drug taken (default: 1)
  --cdose int           count of doses in a regimen (default: 3)
  --inputfile str       specify the inputfile that contains the
                        pharmacokinetic parameters. For deterministic methods
                        please input only one set of pk parameters (default:
                        PK_default.csv)
  --outputfile str      specify the name of output csv file (default: ntm.csv)
  --ifpdf bool          if output the graphic of results in pdf (default:
                        False)
  --ts float            start time point to run the computation, relative to
                        timing of first dose [hr] (default: 0)
  --te float            end time point of the computation, relative to timing
                        of first dose [hr] (default: 240)
  --tspanres float      timespan in the output data [hr] should not be shorter
                        than timesteps (default: 1)
  --ifll bool           if long-lived and latently infected cells are
                        considered (default: False)
  --timesteps number number number
                        time steps for V, T1, T2, [min] Fraction like 1/6 is
                        allowed. For V and T1 the time steps should be same
                        and cannot be shorter than T2 (default: 10 10 1)
  --rate float          percentage of the observed next transition (default:
                        0.999)

```

### ctsm
```
usage: ./prophylaxisEstimator.py ctsm [<options>]

optional arguments:
  -h, --help            show this help message and exit
  --mdose float         mass of DTG dose [mg] (default: 50)
  --adh float           adherence of drug taken (default: 1)
  --cdose int           count of doses in a regimen (default: 3)
  --inputfile str       specify the inputfile that contains the
                        pharmacokinetic parameters. For deterministic methods
                        please input only one set of pk parameters (default:
                        PK_default.csv)
  --outputfile str      specify the name of output csv file (default:
                        ctsm.csv)
  --ifpdf bool          if output the graphic of results in pdf (default:
                        False)
  --ts float            start time point to run the computation, relative to
                        timing of first dose [hr] (default: 0)
  --te float            end time point of the computation, relative to timing
                        of first dose [hr] (default: 240)
  --tspanres float      timespan in the output data [hr] should not be shorter
                        than timesteps (default: 1)
  --ifll bool           if long-lived and latently infected cells are
                        considered (default: False)
  --timesteps number number number
                        time steps for V, T1, T2, [min] Fraction like 1/6 is
                        allowed. For V and T1 the time steps should be same
                        and cannot be shorter than T2 (default: 10 10 1)
```

### pgs

```
usage: ./prophylaxisEstimator.py pgs [<options>]

optional arguments:
  -h, --help        show this help message and exit
  --mdose float     mass of DTG dose [mg] (default: 50)
  --adh float       adherence of drug taken (default: 1)
  --cdose int       count of doses in a regimen (default: 3)
  --inputfile str   specify the inputfile that contains the pharmacokinetic
                    parameters. For deterministic methods please input only
                    one set of pk parameters (default: PK_default.csv)
  --outputfile str  specify the name of output csv file (default: pgs.csv)
  --ifpdf bool      if output the graphic of results in pdf (default: False)
  --ts float        start time point to run the computation, relative to
                    timing of first dose [hr] (default: 0)
  --te float        end time point of the computation, relative to timing of
                    first dose [hr] (default: 240)
  --tspanres float  timespan in the output data [hr] should not be shorter
                    than timesteps (default: 1)
  --ifll bool       if long-lived and latently infected cells are considered
                    (default: False)

```
## Highlights
* compute the prophylactic efficacy with four different methods, including
one stochastic simulation (extrande) and three deterministic methods (ntm, ctsm
and pgs)
* deterministic methods are very efficient and can be applied in many cased,
while extrande can be used to obtain the intermediate state
* return the extinction probabilities and efficacies in csv file
* Uses sane defaults, also provide many possible options

## About
### Running tests
```sh
$ cd test_dir
$ make test
```
#### Clean the generated files after test
```sh
$ make clean
```
