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

(in-package :voxel-octrees)

;; Needed for *voxel*
(deftype dot () '(simple-array single-float))

(declaim (type fixnum *max-dots*))
(defparameter *max-dots* 8
  "Maximum number of dots in leaf")

(declaim (type dot *voxel*))
(defparameter *voxel*
  (make-array 3 :element-type 'single-float :initial-element 1.0)
  "Size of cuboid in set")

(defmacro defvar-unbound (var-name &optional doc)
  `(progn
     (defvar ,var-name)
     (setf (documentation ',var-name 'variable) ,doc)))

(declaim (type (integer 0 #.most-positive-fixnum) *lod*))
(defvar-unbound *lod*
    "May be bound to level of detail. The more LOD is,
     the more accurate result of intersection(s) will
     be calculated. In the other hand the less it is,
     the faster calculations will be. Unbounded *LOD*
     means no level of detail. The least allowed *LOD*
     is zero.")

(defmacro with-lod ((lod) &body body)
  "Stupid wrapper around let.
   Bounds LOD to a certain value."
  `(let ((*lod* ,lod)) ,@body))

(declaim (type (integer 0 #.most-positive-fixnum) *max-depth-local*))
(defparameter *max-depth-local* 2
  "Maximal recursion depth for LOCAL-RAY-TREE-INTERSECTION.
   The more this value is the more farest voxels in space
   are considered `local'. Setting the value too big means
   big penalty in execution time. Try to start with 2.")
