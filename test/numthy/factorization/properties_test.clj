(ns numthy.factorization.properties-test
  (:require [clojure.test :refer :all]
            [numthy.factorization.properties :refer :all]))

(def oeis-1694
  [1, 4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 72, 81, 100, 108, 121, 125, 128, 144,
   169, 196, 200, 216, 225, 243, 256, 288, 289, 324, 343, 361, 392, 400, 432, 441,
   484, 500, 512, 529, 576, 625, 648, 675, 676, 729, 784, 800, 841, 864, 900, 961,
   968, 972, 1000])

(deftest powerful?-test
  (is (= (filter powerful? (range 1 1001))
         oeis-1694)))
