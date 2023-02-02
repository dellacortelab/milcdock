import sys

if len(sys.argv) == 1:
    print("Script to convert coordinates in center/size format to min/max.")
    print("Arg 1: center")
    print("Arg 2: size")
    print("Output: min max")
    sys.exit()

center = float(sys.argv[1])
size = float(sys.argv[2])

minimum = center - size / 2
maximum = center + size / 2

print(f'{minimum:3.3f} {maximum:3.3f}')

