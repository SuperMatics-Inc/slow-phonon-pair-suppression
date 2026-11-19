import Mathlib.Data.Rat.Defs
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

namespace Adiabatic

theorem pairedModes : (2 : ℚ) * (1 / 4) * (1 / 6) = 1 / 12 := by norm_num

theorem wickVertex : (1 / 2 : ℚ) * 4 ^ 2 * (1 / 12) = 2 / 3 := by norm_num

theorem vertexBudget : (2 / 3 : ℚ) * (11 / 4) ^ 8 < 4096 := by norm_num

theorem attachmentAndTemperatures :
    (11 / 4 : ℚ) < 3 ∧ (11 / 4 : ℚ) ^ 2 < 8 ∧ (36 / 5 : ℚ) < 8 ∧
    3 * 16 ^ 2 ≤ (1024 : ℕ) ∧ 8 * 8 ^ 2 ≤ (1024 : ℕ) ∧
    1024 * 4 ^ 3 * 4096 * 64 = 2 ^ 34 ∧
    1024 * (2 ^ 15) ^ 3 * 4096 * 2 * 64 = 2 ^ 74 := by norm_num

theorem tailAndSource :
    4096 ^ 2 / 8 = (2 : ℚ) ^ 21 ∧ 4 * 512 ≤ (4096 : ℚ) ∧
    58 + 6 < (2 : ℚ) ^ 21 - 3 ∧ (1 / 64 : ℚ) + 1 / 64 < 1 / 4 := by norm_num

theorem gaussianAbsorption (y U : ℚ) (hy : 0 ≤ y) (hU : 2048 * U ≤ y) :
    4 * U * (512 * y) ≤ y ^ 2 := by
  nlinarith [mul_nonneg hy (sub_nonneg.mpr hU)]

theorem cutoffLinear (b : ℚ) (hb : 1 ≤ b) :
    1 + 4 * (4097 * b) ≤ 2 ^ 15 * b := by linarith

theorem attachmentLinear (a : ℚ) (ha : 0 ≤ a) :
    4 * (a + 4) ≤ 16 * (1 + 3 * a) := by linarith

theorem attachmentLarge (a : ℚ) (ha : 1 / 2 ≤ a) :
    5 * (4 * (a + 4)) ≤ 36 * (1 + 3 * a) := by linarith

theorem temperatureCube (a b : ℚ) (ha : 0 ≤ a) (hb : 0 ≤ b) :
    a ^ 3 ≤ b ^ 3 ↔ a ≤ b := by
  constructor
  · intro h
    by_contra! hn
    have hsum : 0 < a ^ 2 + a * b + b ^ 2 := by
      have hap : 0 < a := lt_of_le_of_lt hb hn
      nlinarith [mul_pos hap hap, mul_nonneg ha hb, sq_nonneg b]
    nlinarith [mul_pos (sub_pos.mpr hn) hsum]
  · intro h
    have hsum : 0 ≤ b ^ 2 + a * b + a ^ 2 := by
      nlinarith [sq_nonneg a, sq_nonneg b, mul_nonneg ha hb]
    nlinarith [mul_nonneg (sub_nonneg.mpr h) hsum]

theorem quadraticDiscriminant (tau q E : ℚ) :
    (2 * tau - q) ^ 2 - (q ^ 2 + 4 * q * E) = 4 * (tau ^ 2 - q * (E + tau)) := by ring

theorem pairFluxAndJump (l : ℚ) :
    (l / 2) ^ 2 = l ^ 2 / 4 ∧ (1 / 2 : ℚ) * (2 / 2 ^ 2) = 1 / 4 := by
  constructor <;> ring

#print axioms pairedModes
#print axioms wickVertex
#print axioms vertexBudget
#print axioms attachmentAndTemperatures
#print axioms tailAndSource
#print axioms gaussianAbsorption
#print axioms cutoffLinear
#print axioms attachmentLinear
#print axioms attachmentLarge
#print axioms temperatureCube
#print axioms quadraticDiscriminant
#print axioms pairFluxAndJump

end Adiabatic
