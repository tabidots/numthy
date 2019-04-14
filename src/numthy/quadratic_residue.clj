(ns numthy.quadratic-residue
  (:require [clojure.math.numeric-tower :as tower]
            [numthy.modular-arithmetic :as ma]
            [numthy.helpers :as h]
            [numthy.polynomial :refer [degree]]
            [numthy.p-adic :refer [p-adic-order]]))

(defn legendre-symbol
  "Legendre symbol for odd prime moduli. Returns 1 if a is a quadratic residue mod p.
  Note that this uses a probabilistic primality test in order to prevent this
  function from being a bottleneck in other algorithms, when p is sufficiently large."
  [a p]
  (when (and (> p 2) (.isProbablePrime (biginteger p) 5))
    (let [pow (-> (- p 1) (/ 2))  ;; a^((p-1)/ 2) mod p
          sym (ma/mod-pow a pow p)]
      (if (> sym 1) (- sym p)     ;; make all non-residues -1
        sym))))

(defn kronecker-symbol
  "Kronecker symbol, restricted to a positive integer a. Returns 1 if a is
  a quadratic residue mod n."
  ;; https://en.wikipedia.org/wiki/Kronecker_symbol
  [a n]
  (if-let [ls (legendre-symbol a n)] ls  ;; for prime moduli p > 2
    (let [pfs (h/prime-factorization n)] ;; bottleneck for large composite moduli
      (->> (map (fn [p]
                  (if (> p 2) (legendre-symbol a p)
                    (nth [0 1 0 -1 0 -1 0 1] (mod a 8))))
                pfs)
           (reduce *)))))

(defn perfect-power-of-2?
  ;; Maybe this is useful for something
  [n]
  (reduce (fn [a b]
            (cond
              (> a n) (reduced false)
              (= a n) (reduced true)
              :else   (* a b)))
          (repeat 2)))

(defn msqrt
  "Uses the Tonelli-Shanks algorithm to find the integer x s.t. x^2 = a (mod p).
  Returns nil if modulus is composite or if there is no solution."
  ;; https://eli.thegreenplace.net/2009/03/07/computing-modular-square-roots-in-python
  ;; http://www.math.vt.edu/people/brown/doc/sqrts.pdf
  [a p]
  (when-let [ls (legendre-symbol a p)]
    (when (and (ma/coprime? a p) (pos? ls))
            ;; Find a number n that is a quadratic nonresidue mod p
      (let [n (->> (drop 2 (range))
                   (filter #(neg? (legendre-symbol % p)))
                   first)
            ;; Represent (p - 1) as (s * 2^e)
            e (p-adic-order 2 (dec p))
            s (/ (dec p) (reduce *' (repeat e 2)))]
        (loop [x (ma/mod-pow a (-> s (+ 1) (/ 2)) p) ;; 1st guess at sqrt
               b (ma/mod-pow a s p)                  ;; 1st guess at fudge factor
               g (ma/mod-pow n s p)
               r e]
                ;; Find the smallest m < r s.t. b^(2^m) ≡ 1 mod p
          (let [m (->> (range r)
                       (filter #(= 1 (ma/mod-pow b (tower/expt 2 %) p)))
                       first)]
            (if (zero? m) x
              (recur
                (-> (ma/mod-pow g (tower/expt 2 (- r m 1)) p)
                    (* x) (mod p))                         ;; x * g^(2^(r-m-1)) % p
                (-> (ma/mod-pow g (tower/expt 2 (- r m)) p)
                    (* b) (mod p))                         ;; b * g^(2^(r-m)) % p
                (-> (ma/mod-pow g (tower/expt 2 (- r m)) p)
                    (mod p))                               ;; g^(2^(r-m)) % p
                m))))))))

(defn mod-quadratic-zeroes
  [pnml p]
  (when (and (= 2 (degree pnml)) (.isProbablePrime (biginteger p) 5))
    (let [{a 2, b 1, c 0} pnml
          discriminant    (-> (*' b b)
                              (- (* 4 a c)))]
      (when-let [msr-d (msqrt discriminant p)]
        ;; (-b ± sqrt(4ac)) / 2a ==>> (-b ± msqrt(4ac)) * 2a^-1 mod p
        [(-> (- b) (+ msr-d) (* (ma/mod-inverse (* 2 a) p)) (mod p))
         (-> (- b) (- msr-d) (* (ma/mod-inverse (* 2 a) p)) (mod p))]))))
