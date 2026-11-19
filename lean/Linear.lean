import Mathlib.Data.Rat.Defs
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

namespace AdiabaticLinear

theorem frequencyAndWeightedNorm :
    2 * (1 - (1 / 4 : ℚ)) * (1 / 6) = 1 / 4 ∧
    1 + 2 * (1 / 2 : ℚ) * (1 / 4) = 5 / 4 := by norm_num

theorem sumFactor (b : ℚ) :
    3 + 32 * b ≤ 8 * (1 + 4 * b) := by linarith

theorem attachmentAndTemperatures :
    3 * (5 / 4 : ℚ) * 8 < 64 ∧
    1 + 4 * (4097 : ℕ) ≤ 2 ^ 15 ∧
    (64 : ℕ) * 4 * 4096 * 64 = 2 ^ 26 ∧
    (64 : ℕ) * 2 ^ 15 * 4096 * 2 * 64 = 2 ^ 40 ∧
    (64 : ℕ) * 2 ^ 15 * 8 = 2 ^ 24 ∧ 2 ^ 21 - 1 > (30 : ℕ) := by norm_num

theorem inverseTemperature (b t E : ℚ) (hbt : b * t = 1) :
    b ^ 3 * (1 + b * E) * t ^ 4 = E + t := by
  calc
    _ = (b * t) ^ 3 * t + (b * t) ^ 4 * E := by ring
    _ = _ := by rw [hbt]; ring

theorem splitRoot (A t E : ℚ) (hE : A * E ≤ (t ^ 4) / 2)
    (hT : A * t ≤ (t ^ 4) / 2) : A * (E + t) ≤ t ^ 4 := by
  rw [mul_add]
  linarith

theorem stripConstants :
    (1 / 4 : ℚ) ^ 2 / (2 + 1 / 2) = 1 / 40 ∧ (1 / 40 : ℚ) > 1 / 64 ∧
    (1 / 2 : ℚ) ^ 4 = 1 / 16 ∧
    (2 : ℚ) ^ 34 / (128 * 2 ^ 26) ≥ 2 ∧
    (2 : ℚ) ^ 74 / (128 * 2 ^ 40) ≥ 2 ∧
    3 * (3 / 4 : ℚ) / 4 > 1 / 2 := by norm_num

theorem stripFrequencyConstants (L : ℚ) :
    (1 / 4 : ℚ) ^ 2 / (128 * L) = 1 / (2048 * L) ∧
    (1 / 2 : ℚ) ^ 2 / (64 * L * (1 / 2)) = 1 / (128 * L) := by
  constructor
  · norm_num [div_eq_mul_inv, mul_inv_rev]
    ring
  · norm_num [div_eq_mul_inv, mul_inv_rev]
    ring

theorem fillingEndpoint :
    (1 / 32 : ℚ) ^ 2 = 1 / 1024 ∧
    4 * (1 / 32 : ℚ) / 3 < 1 / 4 ∧
    4 * 2 / (4 * 3 * (1 / 32 : ℚ)) < 24 ∧
    (24 : ℕ) < 64 * 2 ^ 26 ∧ (24 : ℕ) < 64 * 2 ^ 40 := by norm_num

#print axioms frequencyAndWeightedNorm
#print axioms sumFactor
#print axioms attachmentAndTemperatures
#print axioms inverseTemperature
#print axioms splitRoot
#print axioms stripConstants
#print axioms stripFrequencyConstants
#print axioms fillingEndpoint
end AdiabaticLinear
