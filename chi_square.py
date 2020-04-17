import sys
sys.path.insert(0,"../")
sys.path.insert(0,"../level_scheme")
sys.path.insert(0,"./")
import trans_strength as ts
test = ts.Vmix(2140, 2332, [0.20, 0.13], [0.055, 0.00], [0.041, .002], error=False, size=100)
print(test)
test1 = ts.Vmix(1947, 2263, [0.30, 0.02], [0.041, 0.003], [0.00, .00], error=False, size=100)
print(test1)
