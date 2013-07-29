;; Copyright (c) 2013, Vasily Postnicov
;; All rights reserved.
;;
;; Redistribution and use in source and binary forms, with or without modification,
;; are permitted provided that the following conditions are met:

;;  Redistributions of source code must retain the above copyright notice, this
;;  list of conditions and the following disclaimer.
;;
;;  Redistributions in binary form must reproduce the above copyright notice, this
;;  list of conditions and the following disclaimer in the documentation and/or
;;  other materials provided with the distribution.
;;
;;THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
;;ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
;;ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
;;(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
;;LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
;;ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

(defpackage voxel-octrees
  (:use #:cl)
  (:export #:make-tree   ; Construct and access tree
           #:node
           #:node-children
           #:node-dots
           #:node-bounding-box
           #:leafp

           #:align-on-voxel ; Geometry functions
           #:calc-abs-metric
           #:calc-sqr-metric
           #:hit-box
           #:hit-plane

           #:ray-tree-intersection ; Tree search
           #:local-ray-tree-intersection

           #:*max-dots* ; Parameters
           #:*voxel*
           #:*max-depth*

           ;; Test values
           #+voxel-octrees-test
           #:*ray-test-search-stats*))
