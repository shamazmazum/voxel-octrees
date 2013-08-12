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

;; (Maybe) faster version

(defun construct-subspace-indices (intersects)
  (let ((res (list (first intersects)))
         (cur-idx (caar intersects)))
    (dolist (inter (cdr intersects))
      (setq cur-idx (logxor cur-idx (ash 1 (car inter))))
      (push (cons cur-idx (cdr inter)) res))
    (reverse res)))

(defun ray-tree-intersection (tree origin dir &optional (depth 0) (path (list tree)))
  "Search intersection of a ray with closest to the ORIGIN
   cuboid in TREE. Returns point in which intersection is occured
   and path to a tree leaf in which intersected voxel is stored."
  (declare (optimize (speed 3))
           (type (integer 0 #.most-positive-fixnum) depth))
  (let* ((bb (node-bounding-box tree))
         (bb-min (car bb))
         (bb-max (cdr bb)))

    (multiple-value-bind (bb-interp bb-inter-coord)
        (and bb (hit-box bb-min bb-max
                         origin dir))
      (cond
        ((not bb-interp) nil)
        
        ((and (boundp '*lod*)
              (= depth *lod*)) (values bb-inter-coord path))
        
        ((leafp tree)
         (let (leaf-intersections
               (dots (node-dots tree)))
           (dotimes (i (length (#+sbcl sb-ext:truly-the
                                #-sbcl the simple-vector dots)))
             (multiple-value-bind (interp coord)
                 (let ((dot (aref dots i)))
                   (hit-box dot (sum-vector (copy-seq (#+sbcl sb-ext:truly-the
                                                       #-sbcl the dot dot)) *voxel*)  ;; Extra allocation
                            origin dir))
               (if interp (push coord leaf-intersections))))
           
           (if leaf-intersections
               (flet ((get-closest (dot1 dot2)
                        (if (< (calc-abs-metric dot1 origin)
                               (calc-abs-metric dot2 origin))
                            dot1 dot2)))
                 (values
                  (reduce #'get-closest leaf-intersections)
                  path)))))
        
        (t
         (let ((center (node-dots tree))
               plane-intersections)
           
           (dotimes (i 3)
             (multiple-value-bind (interp coord) (hit-plane origin dir center i)
               (if (and interp (dot-betweenp coord bb-min bb-max))
                   (push (cons i coord) plane-intersections))))

           (setq plane-intersections (construct-subspace-indices
                                      (cons (cons (get-subspace-idx center bb-inter-coord) bb-inter-coord)
                                            (flet ((get-closest-tagged (dot1 dot2)
                                                     (< (calc-abs-metric (cdr dot1) origin)
                                                        (calc-abs-metric (cdr dot2) origin))))
                                              (sort plane-intersections #'get-closest-tagged)))))

           (dolist (subspace plane-intersections)
             (let ((next-node (nth (car subspace)
                                   (node-children tree))))
               (multiple-value-bind (inter new-path) (ray-tree-intersection next-node
                                                                            (cdr subspace) dir (1+ depth)
                                                                            (cons next-node path))
                 (if inter (return (values inter new-path))))))))))))

;; NB: This DEPTH (which measures degree of locality of two voxel in space) is unrealted
;; to RAY-TREE-INTERSECTION's DEPTH which is related to LOD.

;; TODO: Add example of use, as promised.

;; Excuse me for too long doc string. I feel that I need to make it so
(defun local-ray-tree-intersection (path origin dir &optional (depth 0))
  "Same as RAY-TREE-INTERSECTION but operates not on a tree but on a path
   from a root node to a leaf. When LOCAL-RAY-TREE-INTERSECTION accepts
   path, it performes recursive reverse search for intersection from a leaf
   to a root and stops if intersection is found, PATH is too short or maximum
   recursion level is reached.

   This can be very useful when the PATH is obtained for the first ray with
   origin1 and dir1 and when the same path is used for the second ray which has
   the same origin2=origin1 and dir2 which differs from dir1 only insignificantly.
   (Here comes the prefix LOCAL-)

   The idea for this function comes from the assumption that to `local' rays
   hit two local (in space) voxels which are local in the tree too (voxel octress
   guarantees this). If this is not true (in example in case of set of random voxels),
   be prepared to get penalty in time of execution. See example of use in source code.

   Returns intersection FOUND USING LOCAL PROPERTIES (so if NIL is returned it does not
   mean that there is no intersection at all, try RAY-TREE-INTERSECTION) and a truncated
   PATH to avoid redundant calculations if future calls"
  (declare (optimize (speed 3))
           (type (integer 0 #.most-positive-fixnum) depth))
  (cond
    ((and (/= depth *max-depth-local*)
          (cdr path))
     (let ((inter (ray-tree-intersection (car path) origin dir)))
       (if inter (values inter path)
           (local-ray-tree-intersection (cdr path) origin dir (1+ depth)))))))

;; More precise but slower
;; ?????
#|(defun ray-tree-intersection (tree origin dir)
  "Search intersection of a ray with closest to the ORIGIN
   cuboid in TREE"
  (declare (optimize (speed 3)))

  (cond
    ((leafp tree)
     (let (leaf-intersections
           (dots (node-dots tree)))
       (dotimes (i (length (sb-ext:truly-the simple-vector dots)))
         (multiple-value-bind (interp coord)
             (let ((dot (aref dots i)))
               (hit-box dot (sum-vector (copy-seq (sb-ext:truly-the dot dot)) *voxel*)  ;; Extra allocation
                        origin dir))
           (if interp (push coord leaf-intersections))))

       (if leaf-intersections
           (flet ((get-closest (dot1 dot2)
                    (if (< (calc-abs-metric dot1 origin)
                           (calc-abs-metric dot2 origin))
                        dot1 dot2)))
             (reduce #'get-closest leaf-intersections)))))
        
        (t
         (let ((bb-intersections
                (loop for child in (node-children tree)
                   for bb = (node-bounding-box child)
                   for inter = (if (and bb (hit-box (car bb) (cdr bb) origin dir)) (car bb))
                   when inter collect (cons child inter))))
           
           (setq bb-intersections
                 (flet ((get-closest-tagged (dot1 dot2)
                          (< (calc-abs-metric (cdr dot1) origin)
                             (calc-abs-metric (cdr dot2) origin))))
                   (sort bb-intersections #'get-closest-tagged)))

           (loop for subspace in bb-intersections
              for inter = (ray-tree-intersection (car subspace) origin dir)
              until inter finally (return inter))))))|#
