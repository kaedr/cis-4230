min_number() {
    printf "%sms\n" "$@" | sort -g | head -n1
}

echo "Running 1000x1000"
ONE=$(./GaussianC-VLA.exe ~/1000x1000.dat | cut -d ' ' -f4)
TWO=$(./GaussianC-VLA.exe ~/1000x1000.dat | cut -d ' ' -f4)
THREE=$(./GaussianC-VLA.exe ~/1000x1000.dat | cut -d ' ' -f4)
min_number $ONE $TWO $THREE

echo "Running 1250x1250"
ONE=$(./GaussianC-VLA.exe ~/1250x1250.dat | cut -d ' ' -f4)
TWO=$(./GaussianC-VLA.exe ~/1250x1250.dat | cut -d ' ' -f4)
THREE=$(./GaussianC-VLA.exe ~/1250x1250.dat | cut -d ' ' -f4)
min_number $ONE $TWO $THREE

echo "Running 1500x1500"
ONE=$(./GaussianC-VLA.exe ~/1500x1500.dat | cut -d ' ' -f4)
TWO=$(./GaussianC-VLA.exe ~/1500x1500.dat | cut -d ' ' -f4)
THREE=$(./GaussianC-VLA.exe ~/1500x1500.dat | cut -d ' ' -f4)
min_number $ONE $TWO $THREE

echo "Running 1750x1750"
ONE=$(./GaussianC-VLA.exe ~/1750x1750.dat | cut -d ' ' -f4)
TWO=$(./GaussianC-VLA.exe ~/1750x1750.dat | cut -d ' ' -f4)
THREE=$(./GaussianC-VLA.exe ~/1750x1750.dat | cut -d ' ' -f4)
min_number $ONE $TWO $THREE

echo "Running 2000x2000"
ONE=$(./GaussianC-VLA.exe ~/2000x2000.dat | cut -d ' ' -f4)
TWO=$(./GaussianC-VLA.exe ~/2000x2000.dat | cut -d ' ' -f4)
THREE=$(./GaussianC-VLA.exe ~/2000x2000.dat | cut -d ' ' -f4)
min_number $ONE $TWO $THREE

# echo "Running 5000x5000"
# ./GaussianC-VLA.exe ~/5000x5000.dat
# ./GaussianC-VLA.exe ~/5000x5000.dat
# ./GaussianC-VLA.exe ~/5000x5000.dat

