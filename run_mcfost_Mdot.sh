dump=elias24_01001

rm -rf data_*

m_dot=(0.0 1e-10 3e-10 1e-09 3e-09 1e-08 3e-08 1e-07)
fluffy=(1 0.3 0.1 0.03 0.01)

for mdot in "${m_dot[@]}"
do
    for f in "${fluffy[@]}"
    do
        dir=01001_Mdot_${mdot}_fluffy_${f}

        echo $dir

        /Users/cpinte/mcfost/src/mcfost ../elias24.para -phantom $dump -planet_az -110 -age 1Myr -fluffy $f -mol -freeze-out 20 -photodissociation -photodesorption -Mdot 2 $mdot -Mdot 1 4e-7

        /Users/cpinte/mcfost/src/mcfost ../elias24.para -phantom $dump -planet_az -110 -age 1Myr -fluffy $f -img 3.8  -Mdot 2 $mdot -Mdot 1 4e-7

        /Users/cpinte/mcfost/src/mcfost ../elias24.para -phantom $dump -planet_az -110 -age 1Myr -fluffy $f -img 1300 -Mdot 2 $mdot -Mdot 1 4e-7

        rm -rf $dir
        mkdir -p $dir
        mv data_* $dir
    done
done
