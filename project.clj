(defproject numthy "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :jvm-opts ["-Djdk.attach.allowAttachSelf"]
  :dev-dependencies []
  :dependencies [[org.clojure/clojure "1.10.0"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [org.clojure/math.combinatorics "0.1.4"]
                 [criterium "0.4.4"]
                 [com.clojure-goes-fast/clj-async-profiler "0.3.1"]])
