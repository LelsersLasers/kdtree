-- http://en.wikipedia.org/wiki/K-d_tree

module Data.Trees.KdTree where

import           Data.Maybe

import qualified Data.Foldable   as F
import qualified Data.List       as L
import           Test.QuickCheck

class Point p where
      -- |dimension returns the number of coordinates of a point.
      dimension :: p -> Int

      -- |coord gets the k'th coordinate, starting from 0.
      coord :: Int -> p -> Double

      -- |dist2 returns the squared distance between two points.
      -- dist2 :: p -> p -> Double
      -- dist2 a b = sum . map diff2 $ [0..dimension a - 1]
      --   where diff2 i = (coord i a - coord i b)^2

-- |compareDistance p a b  compares the distances of a and b to p.
-- compareDistance :: (Point p) => p -> p -> p -> Ordering
-- compareDistance p a b = dist2 p a `compare` dist2 p b

data Point3d = Point3d { p3x :: Double, p3y :: Double, p3z :: Double }
    deriving (Eq, Ord, Show)

instance Point Point3d where
    dimension _ = 3

    coord 0 p = p3x p
    coord 1 p = p3y p
    coord 2 p = p3z p


data KdTree point =
    KdNode {
        kdLeft  :: KdTree point,
        kdPoint :: point,
        kdRight :: KdTree point,
        kdAxis  :: Int
    }
    | KdEmpty deriving (Eq, Ord, Show)

instance Functor KdTree where
    fmap _ KdEmpty             = KdEmpty
    fmap f (KdNode l x r axis) = KdNode (fmap f l) (f x) (fmap f r) axis

instance F.Foldable KdTree where
    foldr f init KdEmpty = init
    foldr f init (KdNode l x r _) = F.foldr f init3 l
        where init3 = f x init2
              init2 = F.foldr f init r

fromList :: Point p => [p] -> KdTree p
fromList points = fromListWithDepth points 0

-- |fromListWithDepth selects an axis based on depth so that the axis cycles
-- through all valid values.
fromListWithDepth :: Point p => [p] -> Int -> KdTree p
fromListWithDepth [] _ = KdEmpty
fromListWithDepth points depth = node
    where axis = depth `mod` dimension (head points)

          -- Sort point list and choose median as pivot element
          sortedPoints =
              L.sortBy (\a b -> coord axis a `compare` coord axis b) points
          medianIndex = length sortedPoints `div` 2
          medianCoordinate = coord axis (sortedPoints !! medianIndex)

          leftPoints = filter (\p -> coord axis p < medianCoordinate) sortedPoints
          trueMedianIndex = length leftPoints
          rightPoints = drop (trueMedianIndex+1) sortedPoints

          -- Create node and construct subtrees
          node = KdNode { kdLeft  = fromListWithDepth leftPoints (depth+1),
                          kdPoint = sortedPoints !! trueMedianIndex,
                          kdRight = fromListWithDepth rightPoints (depth+1),
                          kdAxis  = axis }

toList :: KdTree p -> [p]
toList = F.foldr (:) []

-- |subtrees t returns a list containing t and all its subtrees, including the
-- empty leaf nodes.
subtrees :: KdTree p -> [KdTree p]
subtrees KdEmpty               = [KdEmpty]
subtrees t@(KdNode l x r axis) = subtrees l ++ [t] ++ subtrees r

-- |nearestNeighbor tree p returns the nearest neighbor of p in tree using Euclidean distance
nearestNeighbor :: Point p => KdTree p -> p -> Maybe p
nearestNeighbor = nearestNeighborBy dist
  where
    dist a b = sum [ (coord i a - coord i b) ^ 2 | i <- [0 .. dimension a - 1] ]

-- nearestNeighbor KdEmpty probe = Nothing
-- nearestNeighbor (KdNode KdEmpty p KdEmpty _) probe = Just p
-- nearestNeighbor (KdNode l pivot r axis) probe =
--     if xProbe < xPivot then findNearest l r else findNearest r l
--     where xProbe = coord axis probe
--           xPivot = coord axis pivot
--           findNearest tree1 tree2 =
--                 let candidate1 = case nearestNeighbor tree1 probe of
--                                    Nothing   -> pivot
--                                    Just best -> L.minimumBy (compareDistance probe) [best, pivot]
--                     sphereIntersectsPlane = (xProbe - xPivot)^2 <= dist2 probe candidate1
--                     candidates2 = if sphereIntersectsPlane
--                                     then candidate1 : maybeToList (nearestNeighbor tree2 probe)
--                                     else [candidate1] in
--                 Just . L.minimumBy (compareDistance probe) $ candidates2


-- |nearestNeighbor tree p returns the nearest neighbor of p in tree using the distance function
nearestNeighborBy :: (Point p) => (p -> p -> Double) -> KdTree p -> p -> Maybe p
nearestNeighborBy _ KdEmpty _ = Nothing
nearestNeighborBy _ (KdNode KdEmpty p KdEmpty _) _ = Just p
nearestNeighborBy dist (KdNode l pivot r axis) probe =
    if xProbe < xPivot then search l r else search r l
  where
    xProbe = coord axis probe
    xPivot = coord axis pivot
    search tree1 tree2 =
        let best1 = fromMaybe pivot (nearestNeighborBy dist tree1 probe)
            bestSoFar = if dist probe best1 < dist probe pivot then best1 else pivot
            shouldCheckOtherSide = (xProbe - xPivot)^2 <= dist probe bestSoFar
            bests = if shouldCheckOtherSide
                      then bestSoFar : maybeToList (nearestNeighborBy dist tree2 probe)
                      else [bestSoFar]
        in Just $ L.minimumBy (\a b -> compare (dist probe a) (dist probe b)) bests


-- |nearNeighbors tree p returns all neighbors within distance r from p in tree.
nearNeighbors :: Point p => (p -> p -> Double) -> KdTree p -> Double -> p -> [p]
nearNeighbors dist KdEmpty radius probe                      = []
nearNeighbors dist (KdNode KdEmpty p KdEmpty _) radius probe = [p | dist p probe <= radius^2]
nearNeighbors dist (KdNode l p r axis) radius probe          =
    if xProbe <= xp
      then let nearest = maybePivot ++ nearNeighbors dist l radius probe
           in if xProbe + abs radius > xp
                then nearNeighbors dist r radius probe ++ nearest
                else nearest
      else let nearest = maybePivot ++ nearNeighbors dist r radius probe
           in if xProbe - abs radius < xp
                then nearNeighbors dist l radius probe ++ nearest
                else nearest
  where xProbe     = coord axis probe
        xp         = coord axis p
        maybePivot = [p | dist probe p <= radius^2]

-- |isValid tells whether the K-D tree property holds for a given tree.
-- Specifically, it tests that all points in the left subtree lie to the left
-- of the plane, p is on the plane, and all points in the right subtree lie to
-- the right.
isValid :: Point p => KdTree p -> Bool
isValid KdEmpty = True
isValid (KdNode l p r axis) = leftIsGood && rightIsGood
    where x = coord axis p
          leftIsGood = all ((<= x) . coord axis) (toList l)
          rightIsGood = all ((>= x) . coord axis) (toList r)

-- |allSubtreesAreValid tells whether the K-D tree property holds for the given
-- tree and all subtrees.
allSubtreesAreValid :: Point p => KdTree p -> Bool
allSubtreesAreValid = all isValid . subtrees

-- |kNearestNeighbors tree k p returns the k closest points to p within tree.
kNearestNeighbors :: (Eq p, Point p) => (p -> p -> Double) -> KdTree p -> Int -> p -> [p]
kNearestNeighbors _ KdEmpty _ _ = []
kNearestNeighbors _ _ k _ | k <= 0 = []
kNearestNeighbors dist tree k probe = nearest : kNearestNeighbors dist tree' (k-1) probe
    where nearest = fromJust $ nearestNeighborBy dist tree probe
          tree' = tree `remove` nearest

-- |remove t p removes the point p from t.
remove :: (Eq p, Point p) => KdTree p -> p -> KdTree p
remove KdEmpty _ = KdEmpty
remove (KdNode l p r axis) pKill
  | p == pKill = fromListWithDepth (toList l ++ toList r) axis
  | coord axis pKill <= coord axis p = KdNode (remove l pKill) p r axis
  | otherwise = KdNode l p (remove r pKill) axis

instance Arbitrary Point3d where
    arbitrary = do
        x <- arbitrary
        y <- arbitrary
        z <- arbitrary
        return (Point3d x y z)

