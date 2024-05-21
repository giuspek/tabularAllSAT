#!/bin/sh

solver="solver"
dualiza="../dualiza/dualiza"

if [ ! -x $solver ]
then
  echo "test.sh: error: could not find tabularasat"
  exit 1
fi

if [ ! -x "$dualiza" ]
then
  echo "test.sh: error: could not find dualiza"
  exit 1
fi

run () {
  cnf=$1

  # Run solvers.
  echo "Solving $cnf ..."
  cmd="./$dualiza $cnf -c"
  # sharpsat_mc=$(eval $cmd | sed -e '1,/# solutions/d' -e '/# END/,$d')
  sharpsat_mc=$(eval $cmd | tail -1)

  cmd1="./$solver $cnf -q -c"
  gcai19_mc=$(eval $cmd1 | tail -1)

  status=$?
  if [ ! $gcai19_mc = $sharpsat_mc ]
  then
    echo "test.sh: error: dualiza and mc_tabularasat(NON_CHRONO_PRIME) do not match on $cnf"
    echo "dualiza: $sharpsat_mc"
    echo "MC_TABULARASAT: $gcai19_mc"
  fi

}

for i in ../test/basic/*.cnf
do
  run $i
done

exit 0
