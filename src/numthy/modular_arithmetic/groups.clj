(ns numthy.modular-arithmetic.groups
  (:require [clojure.math.numeric-tower :refer [expt lcm]]
            [numthy.factorization.core :refer [factorize phi]]
            [numthy.helpers :refer [coprime?]]
            [numthy.modular-arithmetic.utils :refer [mod-mul mod-pow odd-prime-power?]]
            [numthy.perfect-powers :refer [perfect-powers]]
            [numthy.primes.is-prime :refer [prime?]]))

(defn complete-residue-system?
  "A complete residue system mod m is a set containing precisely one
  representative of each residue class (congruence class) mod m.
  i.e., {-13 4 17 18} mod 4 = {3 0 1 2} ∴ true"
  [coll m]
  (= (set (range m)) (set (map #(mod % m) coll))))

(defn generators-additive ;; https://www.doc.ic.ac.uk/~mrh/330tutor/ch06.html
  "Generators of (ℤ/nℤ,+), the additive group of integers mod n—
  i.e., the integers less than n that generate all integers < n
  when added repeatedly to themselves, excluding 0 (universal non-generator) and
  1 (universal generator)."
  [n]
  (let [group    (set (range n))
        generate (fn [seed] (set (map #(mod-mul seed % n) (range n))))]
    (->> (drop 2 (range n))
         (filter #(= group (generate %))))))

(defn multiplicative-group
  "(ℤ/nℤ,×), the multiplicative group of integers mod n,
  or all the integers < n that have a multiplicative inverse mod n.
  Effectively, it is all of the positive integers less than n that are coprime to n."
  [n]
  (if (prime? n) (range 1 n)
    (filter #(coprime? n %) (range 1 n))))

(comment
 "Multiplicative order also enables an alternative test for primitive root-ness.
  Works for prime and composite moduli, but will be impractically slow for large moduli."
 (= (multiplicative-order x m) (h/phi p)))

(defn cyclic?
  "Returns true if (ℤ/nℤ,×) ≅ C_Φ(n), i.e., the multiplicative group of integers mod n
  is cyclic. That is, the group has at least one primitive root. Returns nil otherwise."
  ;; This uses Gauss's proof that (ℤ/nℤ,×) is cyclic iff n = 1, 2, 4, p^k, or 2p^k, p prime and k positive.
  [n]
  (or (contains? #{1 2 4} n)
      (odd-prime-power? n)
      (odd-prime-power? (* n 1/2))))

(comment
 "Naive methods to determine whether a group is cyclic"
 (when (.isProbablePrime (biginteger n) 5)
   (some #(primitive-root? % n) (range 1 n)))
 (some? (primitive-roots n)))

(defn carmichael
  "The Carmichael function of n gives the exponent of (ℤ/nℤ,×), i.e. the smallest positive
  integer k such that a^k ≡ 1 (mod n) for all a coprime to n."
  [n]
  (if (= n 1) 1
    (letfn [(prime-power-carm [prime power]
              ;; for a power of an odd prime and for 2 and 4, λ(n) is equal to
              ;; the Euler totient φ(n); for powers of 2 greater than 4 it is
              ;; equal to half of the Euler totient.
              (cond
                (and (odd? prime)
                     (= 1 power)) (dec prime) ;; speed up calculation
                (odd? prime)      (phi (expt prime power))
                (> power 2)       (/ (phi (expt prime power)) 2)
                :else             (phi (expt prime power))))]
      ;; λ(n) is the least common multiple of the λ of each of its prime power factors
      (reduce-kv (fn [res prime power]
                   (lcm res (prime-power-carm prime power)))
                 1 (factorize n)))))

(comment
 "Alternative check for cyclic-ness"
 (= (phi n) (carmichael n)))

;; TODO: https://en.wikipedia.org/wiki/Carmichael_number
;; TODO: https://en.wikipedia.org/wiki/Euler_pseudoprime
;; https://en.wikipedia.org/wiki/Fermat%27s_little_theorem

;; https://en.wikipedia.org/wiki/Congruent_number

;; TODO: Dirichlet character https://en.wikipedia.org/wiki/Dirichlet_character
