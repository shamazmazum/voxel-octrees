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

(declaim (inline leafp))
(defun leafp (node)
  "Is node a leaf?"
  (null (node-children node)))

(defstruct (node (:print-function (lambda (struct stream depth)
                                    (declare (ignore depth))
                                    (print-unreadable-object (struct stream
                                                              :type t
                                                              :identity t)
                                      (format stream "LEAFP:~A" (leafp struct))))))
  children
  bounding-box
  dots)

(declaim (ftype (function (dot dot &optional dot) dot) sum-vector))
(defun sum-vector (dot1 dot2 &optional (res (copy-seq dot1)))
  "Returns sum of two dots. If RES is supplied, destructively
   modifies it rather then creating new dot."
  (declare (type dot dot1 dot2 res)
           (optimize (speed 3)))
  (map-into res #'+ dot1 dot2))

;; Commented ones are more readable but slower examples

#|(defun calc-bounding-box (dots)
  (cons (reduce #'(lambda (x y) (map 'vector #'min x y)) dots)
        (sum-vector
         (reduce #'(lambda (x y) (map 'vector #'max x y)) dots)
         *voxel*)))|#

(defun calc-bounding-box (dots)
  "Calc minimal cuboid hull for set of cuboids"
  (declare (type simple-vector dots)
           (optimize (speed 3)))
  
  (let* ((first-dot (aref dots 0))
         (min (copy-seq (the dot first-dot)))
         (max (copy-seq (the dot first-dot))))
    (loop for x across dots do
         (map-into min #'min min (the dot x))
         (map-into max #'max max (the dot x)))
    (cons min (sum-vector max *voxel* max))))


#|(defun calc-avg (dots)
  (map 'vector
       #'(lambda (x) (/ x (length dots)))
       (reduce #'(lambda (x y) (map 'vector #'+ x y)) dots)))|#

(defun calc-avg (dots)
  "Calculate dot in center of set"
  (declare (type simple-vector dots)
           (optimize (speed 3)))
  (let* ((sum-start (make-array 3
                               :element-type 'single-float
                               :initial-element 0.0))
         (sum (reduce #'(lambda (d1 d2)
                          (sum-vector d1 d2 d1))
                      dots :initial-value sum-start)))
    
    (declare (type dot sum))
    (map-into sum #'(lambda (x)
                      (declare (type single-float x))
                      (/ x (length dots)))
              sum)))

(declaim (ftype (function (dot dot) fixnum) get-subspace-idx))
(defun get-subspace-idx (dot1 dot2)
  "Placement of dot is encoded by subspace index
   In 3D there are 2^3 different possibilities
   of placement around the center. This function
   calculates which possibility is the case."
  (declare (optimize (speed 3)))

  (let ((flags (make-array 3 :element-type '(unsigned-byte 8))))
    (declare (type dot dot1 dot2)
             (type (simple-array (unsigned-byte 8)) flags)
             (dynamic-extent flags))

    (flet ((get-flag (c1 c2 i)
             (declare (type fixnum i)
                      (type single-float c1 c2))
             (ash (if (> c1 c2) 1 0) i)))
      
      (reduce #'logior (map-into flags #'get-flag dot1 dot2 '(0 1 2))))))
  
(defun align-on-voxel (dot)
  "Destructevily aligns on voxel"
  (declare (optimize (speed 3))
           (type dot dot))

  (flet ((align-coord (x alignment)
           (declare (type single-float x alignment))
           (* alignment (fceiling x alignment))))
    (map-into dot #'align-coord dot *voxel*)))

;; Each axis-aligned voxel is represented by its minimal coordinate
;; and its sides stored in *VOXEL*.
;; (For purposes this library was designed for, all voxels has
;; the same length).

(defun make-tree (dots)
  "Build tree based on set of axis-aligned cuboids with min coordinates
   stored in dots"
  (declare (type simple-vector dots)
           (optimize (speed 3)))
  (let ((node (make-node :bounding-box (and (/= (length dots) 0)
                                            (calc-bounding-box dots)))))
      
    (if (< (length dots) *max-dots*) (setf (node-dots node) dots)
        (let ((avg (calc-avg dots)))
          (setf (node-dots node) (align-on-voxel avg)
                (node-children node)
                (loop for i fixnum below 8 collect
                     (flet ((current-subspace (dot)
                              (= i (get-subspace-idx avg dot))))
                       (make-tree (remove-if-not #'current-subspace dots)))))))
    node))

(defun inaccurate-tree-depth (tree &optional (depth 0))
  "Returns inaccurate depth of tree. Leaf node has zero
   depth."
  (declare (optimize (speed 3))
           (type (integer 0 #.most-positive-fixnum) depth))
  (if (leafp tree) depth
      (inaccurate-tree-depth (nth
                              (logand #x07 depth)
                              (node-children tree))
                             (1+ depth))))

(defun voxels-in-tree (tree)
  "Returns number of voxels in tree"
  (if (leafp tree) (length (node-dots tree))
      (reduce #'+ (mapcar #'voxels-in-tree (node-children tree)))))
