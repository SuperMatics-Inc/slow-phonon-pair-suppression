import Mathlib.Data.Rat.Defs
import Mathlib.Tactic.NormNum

namespace AdiabaticWitness

theorem spectralSchur :
    (1 : ℚ) + (5/12 + (3/4+4/9)/5)/(5/4) < (2*(157/50)*(1/5))^2 ∧
    (-17/27 : ℚ) + 25/16 + 2/5 < (2*(157/50)*(1/5))^2 ∧
    (65/144 : ℚ) + 17/48 + 2/3 + 1/20 < (2*(157/50)*(1/5))^2 := by norm_num

theorem spectralTrial :
    ((1 : ℚ) + 2*(2/3)/3*(5/4) + 2*(1/4)/5*(13/36) +
      ((2/3)/3)^2*(-17/9) + 2*((2/3)/3)*((1/4)/5)*(17/16) +
      ((1/4)/5)^2*(-123/50)) / (1 + (2/3)^2/3 + (1/4)^2/5) >
    (2*(22/7)*(9/50))^2 := by norm_num

theorem dyadicProducts :
    (9/50 : ℚ)^2 * (99/1000) * (24000 - 8) = 24048981/312500 ∧
    (1/5 : ℚ)^2 * (101/900) * (25200 + 8) = 636502/5625 ∧
    (64 : ℚ) < 24048981/312500 ∧ (636502/5625 : ℚ) < 128 := by norm_num

theorem ivtWitness : (96 : ℚ) * (3/4) / 512 = (3/8)^2 := by norm_num

theorem logAndDOS :
    (2/3 : ℚ) + (2/3)/27/(1-1/9) = 25/36 ∧ (25/36 : ℚ) < 7/10 ∧
    (99/1000 : ℚ) < 1 / ((22/7)^2 * (1+1/1024)) ∧
    ((1 : ℚ) + (1/1024)^2)/9 < 101/900 := by norm_num

#print axioms spectralSchur
#print axioms spectralTrial
#print axioms dyadicProducts
#print axioms ivtWitness
#print axioms logAndDOS
end AdiabaticWitness
