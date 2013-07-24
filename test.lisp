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

(defpackage voxel-octrees-test
  (:use #:cl #:voxel-octrees)
  (:export #:verify-ray-tree1
           #:verify-ray-tree2

           #:check-tree))

(in-package :voxel-octrees-test)

(defun construct-set (n &optional (side-size 100000.0))
  (declare (type single-float side-size))
  (let ((array (make-array (list n))))
    (dotimes (i n)
      (setf (aref array i)
            (align-on-voxel
             (make-array 3
                         :element-type 'single-float
                         :initial-contents (list (- (random side-size)
                                                    (/ side-size 2))
                                                 
                                                 (+ 300.0 (random side-size))

                                                 (- (random side-size)
                                                    (/ side-size 2)))))))
    array))

(defun naive-search (set origin dir)
  (let (intersections)
    (dotimes (i (length set))
      (multiple-value-bind (interp coord)
          (let ((dot (aref set i)))
            (hit-box dot
                     (voxel-octrees::sum-vector (copy-seq dot) *voxel*)
                     origin dir))
        (if interp (push coord intersections))))

    (flet ((get-closest (dot1 dot2)
             (if (< (calc-sqr-metric origin dot1)
                    (calc-sqr-metric origin dot2))
                 dot1 dot2)))
      (if intersections
          (reduce #'get-closest intersections)))))

(defun verify-ray-tree1 (n-rays)
  "Construct large set of dots and scan through it
   both with naive O(n) search and with fast ray-tree search.
   Rays of number N-RAYS are emitted from a single start point"
  (let* ((origin (make-array 3
                             :element-type 'single-float
                             :initial-element 0.0))
         (hit-count 0)
         (side-size 2000.0)
         (set (construct-set 500000 side-size))
         (tree (make-tree set))
         (allowed-distance (calc-sqr-metric origin *voxel*)))
    
    (dotimes (i n-rays)
      (let* ((dir (make-array 3
                              :element-type 'single-float
                              :initial-contents (list
                                                 (- (random side-size)
                                                    (/ side-size 2)
                                                    10.0)
                                                 280.0
                                                 (- (random side-size)
                                                    (/ side-size 2)
                                                    10.0))))
             (inter1 (naive-search set origin dir))
             (inter2 (ray-tree-intersection tree origin dir)))

        (if (or
             (or (and inter1 (not inter2))
                 (and inter2 (not inter1)))
             
             (and inter1 inter2
                  (> (calc-sqr-metric inter1 inter2)
                     allowed-distance)))
            
            (return-from verify-ray-tree1
              (progn
                (format t "RAY-TREE-INTERSECTION mismatch: ~A against ~A.~%Was ~D successful hits before fail~%"
                        inter1 inter2 hit-count)
                nil)))
        (if inter1 (incf hit-count)))
      (format t "Ray ~D/~D~%" i n-rays))
    
    (format t "Test succeeds. ~D hits total~%" hit-count)))

(defun verify-ray-tree2 (n-rays)
  "Construct large set of dots and scan through it
   both with naive O(n) search and with fast ray-tree search.
   Rays of number N-RAYS are emitted from different start points, but
   have common direction"
  (let* ((dir (make-array 3
                             :element-type 'single-float
                             :initial-contents '(10.0 10.0 0.0)))
         (hit-count 0)
         (side-size 2000.0)
         (set (construct-set 500000 side-size))
         (tree (make-tree set))
         (allowed-distance (calc-sqr-metric (make-array 3
                             :element-type 'single-float
                             :initial-element 0.0) *voxel*)))
    
    (dotimes (i n-rays)
      (let* ((origin (make-array 3
                                 :element-type 'single-float
                                 :initial-contents (list
                                                    (- (random side-size)
                                                       (/ side-size 2)
                                                       300.0)
                                                    
                                                    0.0
                                                    
                                                    (- (random side-size)
                                                       (/ side-size 2)
                                                       10.0))))
             (inter1 (naive-search set origin dir))
             (inter2 (ray-tree-intersection tree origin dir)))

        (if (or
             (or (and inter1 (not inter2))
                 (and inter2 (not inter1)))
             
             (and inter1 inter2
                  (> (calc-sqr-metric inter1 inter2)
                     allowed-distance)))
            
            (return-from verify-ray-tree2
              (progn
                (format t "RAY-TREE-INTERSECTION mismatch: ~A against ~A.~%Was ~D successful hits before fail~%"
                        inter1 inter2 hit-count)
                nil)))
        (if inter1 (incf hit-count)))
      (format t "Ray ~D/~D~%" i n-rays))
    
    (format t "Test succeeds. ~D hits total~%" hit-count)))

(defun check-tree (tree)
  "Check if tree is constructed correctly"
  (let* ((bb (node-bounding-box tree))
         (bb-min (car bb))
         (bb-max (cdr bb)))
    (cond
      ((leafp tree)
       (if bb
           ;; Check is all voxels are contained in bounding-box
           (every #'(lambda (dot)
                      (and (voxel-octrees::dot-betweenp dot bb-min bb-max)
                           (voxel-octrees::dot-betweenp (voxel-octrees::sum-vector (copy-seq dot) *voxel*) bb-min bb-max)))
                  (node-dots tree))
           t))

      (t
       (let ((center (node-dots tree)))
         (every #'(lambda (x) x) ;; FIXME: This is ugly
                (loop for child in (node-children tree)
                   for i from 0 by 1 collect
                     (let ((child-bb (node-bounding-box child)))
                       (and
                        (if child-bb
                            (and
                             ;; Check if bounding box lays in correct subspace
                             (= i (voxel-octrees::get-subspace-idx center (car child-bb)))
                             ;; FIXME: add something about (cdr child-bb)
                        
                             ;; Check if child's bounding box is contained by parent's
                             (voxel-octrees::dot-betweenp (car child-bb) bb-min bb-max)
                             (voxel-octrees::dot-betweenp (cdr child-bb) bb-min bb-max))
                            t)

                        ;; Recursively check children
                        (check-tree child))))))))))
