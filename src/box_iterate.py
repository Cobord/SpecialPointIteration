"""
iterate within a prod_i [a_i,b_i]
"""

from enum import Enum, auto
from typing import Callable, List, Optional, Tuple
from math import sin, pi as PI, asin
# pylint:disable=import-error
from adic_rational import AdicRational
from parameter_iterate import MultiParameterIterator

class BoxTransformation(Enum):
    """
    how to specify the transfomations from
    the standard [-1,1]^d to the provided box
    """
    NOTRANSFORMATION = auto()
    COORDINATEWISELINEAR = auto()
    COORDINATEWISESIN = auto()
    COORDINATEWISEASIN = auto()
    COORDINATEWISEGIVEN = auto()
    CUSTOM = auto()

class BoxIterator:
    """
    iterate through special points in prod_i [a_i,b_i]
    transformed from the standard one which goes through
    special points in [-1,1]^d
    as described in parameter_iterate.py
    """
    #pylint:disable=too-many-locals,too-many-branches
    def __init__(self, *, denominator_localization : int = 2,
                 num_dimensions : Optional[int] = None,
                 bounding_box : Optional[List[Tuple[float,float]]]=None,
                 expected_per_denominator : Optional[Callable[[int],Optional[int]]] = None,
                 hard_denominator_power_cut : Optional[int] = None,
                 transformation : BoxTransformation = BoxTransformation.COORDINATEWISELINEAR,
                 coordinatewise_given_func : Optional[Callable[[int,float],float]]=None,
                 custom_func : Optional[Callable[[List[float]],List[float]]]=None,
                 bdry_checking_eps : float = 1e-3
                 ):
        bounding_box = self.__check_bounding_box(bounding_box,num_dimensions)
        self.__pre_transformation_iter = MultiParameterIterator(
                denominator_localization,
                self.__num_dimensions,
                expected_per_denominator,
                hard_denominator_power_cut)
        self.__how_many_coarsest(expected_per_denominator,hard_denominator_power_cut)
        if transformation == BoxTransformation.NOTRANSFORMATION:
            for (left_pt,right_pt) in bounding_box:
                if abs(left_pt+1.0)>bdry_checking_eps or abs(right_pt-1.0)>bdry_checking_eps:
                    raise ValueError("When no transform, bounding box should be [-1,1]^d")
            self.post_transformation : Callable[[Tuple[AdicRational,...]],Tuple[float,...]] = \
                lambda x: tuple(map(lambda z: z.to_float(),x))
        elif transformation == BoxTransformation.COORDINATEWISELINEAR:
            def linear_transform(standard_pt : Tuple[AdicRational,...]) -> Tuple[float,...]:
                new_point = [0.0]*len(standard_pt)
                for (idx,coord) in enumerate(standard_pt):
                    left_bdry, right_bdry = bounding_box[idx]
                    new_point[idx] = (coord.to_float() + 1.0)/2.0*(right_bdry-left_bdry) + left_bdry
                return tuple(new_point)
            self.post_transformation = linear_transform
        elif transformation == BoxTransformation.COORDINATEWISESIN:
            def sin_transform(standard_pt : Tuple[AdicRational,...]) -> Tuple[float,...]:
                new_point = [0.0]*len(standard_pt)
                for (idx,coord) in enumerate(standard_pt):
                    left_bdry, right_bdry = bounding_box[idx]
                    new_point[idx] = (sin(PI/2.0*coord.to_float()) + 1.0)/2.0*\
                        (right_bdry-left_bdry) + left_bdry
                return tuple(new_point)
            self.post_transformation = sin_transform
        elif transformation == BoxTransformation.COORDINATEWISEASIN:
            def asin_transform(standard_pt : Tuple[AdicRational,...]) -> Tuple[float,...]:
                new_point = [0.0]*len(standard_pt)
                for (idx,coord) in enumerate(standard_pt):
                    left_bdry, right_bdry = bounding_box[idx]
                    new_point[idx] = (asin(coord.to_float())*2.0/PI + 1.0)/2.0*\
                        (right_bdry-left_bdry) + left_bdry
                return tuple(new_point)
            self.post_transformation = asin_transform
        elif transformation == BoxTransformation.COORDINATEWISEGIVEN:
            if coordinatewise_given_func is None:
                raise ValueError("A coordinatewise transformation must be given")
            # only check the endpoints
            # otherwise
            # assuming the caller actually gave a [-1,1] simeq [a_i,b_i]
            BoxIterator.__check_boundaries(bounding_box,coordinatewise_given_func,bdry_checking_eps)
            def given_transform(standard_pt : Tuple[AdicRational,...]) -> Tuple[float,...]:
                new_point = [0.0]*len(standard_pt)
                for (axis_num,coord) in enumerate(standard_pt):
                    new_point[axis_num] = coordinatewise_given_func(axis_num,coord.to_float())
                return tuple(new_point)
            self.post_transformation = given_transform
        elif transformation == BoxTransformation.CUSTOM:
            if custom_func is None:
                raise ValueError("A custom transformation must be given")
            # does not check anything about custom_func
            # assuming the caller actually gave a [-1,1]^d simeq prod_i [a_i,b_i]
            def given_custom_transform(standard_pt : Tuple[AdicRational,...]) -> Tuple[float,...]:
                standard_pt_floated = list(map(lambda x : x.to_float(),standard_pt))
                new_point = custom_func(standard_pt_floated)
                return tuple(new_point)
            self.post_transformation = given_custom_transform
        else:
            assert False, "A transformation must be one of the enum"

    def my_bounding_box(self) -> List[Tuple[float,float]]:
        """
        give a copy of the bounding box
        """
        return list(self.__bounding_box)

    def __iter__(self):
        """
        itself
        """
        return self

    def __next__(self) -> Tuple[float,...]:
        """
        the next point as a tuple of num_dimensions floats
        """
        standard_pt = self.__pre_transformation_iter.__next__()
        return self.post_transformation(standard_pt)

    def __check_bounding_box(self,bounding_box : Optional[List[Tuple[float,float]]],
                             num_dimensions : Optional[int]) -> None:
        if bounding_box is None and num_dimensions is None:
            raise ValueError(
                "At least the number of dimensions or the bounding box must be explicit")
        if bounding_box is None:
            bounding_box = [(0.0,1.0) for _ in range(num_dimensions)]
        if len(bounding_box) < 1:
            raise ValueError("There must be at least one coordinate axis")
        for (left_pt,right_pt) in bounding_box:
            if right_pt<=left_pt:
                raise ValueError("Each axis of the box must be an interval")
        self.__num_dimensions = len(bounding_box)
        self.__bounding_box = bounding_box
        return bounding_box

    @staticmethod
    def __check_boundaries(bounding_box : List[Tuple[float,float]],
                           coordinatewise_given_func : Callable[[int,float],float],
                           bdry_checking_eps : float) -> None:
        """
        helper for coordinatewise given function
        making sure it maps the -1,1 for each axis
        to the respective a_i and b_i
        """
        for axis_num,(left_bdry,right_bdry) in enumerate(bounding_box):
            neg_one_goes = coordinatewise_given_func(axis_num,-1.0)
            one_goes = coordinatewise_given_func(axis_num,1.0)
            seen_left_endpoint, seen_right_endpoint = \
                min(neg_one_goes,one_goes), max(neg_one_goes,one_goes)
            if abs(seen_left_endpoint - left_bdry)>bdry_checking_eps or \
                abs(seen_right_endpoint - right_bdry)>bdry_checking_eps:
                #pylint:disable=line-too-long
                raise ValueError(
                    f"On axis {axis_num} the transformation should take -1,1 to {left_bdry},{right_bdry} or vice-versa")

    def __how_many_coarsest(self,expected_per_denominator,hard_denominator_power_cut):
        """
        helper for number of points where
        the coordinates pre-transformation
        are all -1,0,1
        """
        if expected_per_denominator is not None:
            self.__num_with_denom_power_0 = expected_per_denominator(0)
        else:
            self.__num_with_denom_power_0 = None
        if self.__num_with_denom_power_0 is None:
            self.__num_with_denom_power_0 = 3**self.__num_dimensions
        if hard_denominator_power_cut is not None and hard_denominator_power_cut<0:
            self.__num_with_denom_power_0 = 0

if __name__ == "__main__":
    np = BoxIterator(denominator_localization=2,bounding_box=[(0,1),(-1,1)])
    # pylint:disable = invalid-name
    index = 0
    for cur in np:
        print(cur)
        index += 1
        if index>20:
            break
    print("Now with sin")
    np = BoxIterator(denominator_localization=2,bounding_box=[(0,1),(-1,1)],
                     transformation=BoxTransformation.COORDINATEWISESIN)
    # pylint:disable = invalid-name
    index = 0
    for cur in np:
        print(cur)
        index += 1
        if index>20:
            break
    print("Now with asin")
    np2 = BoxIterator(denominator_localization=2,
                      num_dimensions=1,
                     transformation=BoxTransformation.COORDINATEWISEASIN)
    # pylint:disable = invalid-name
    index = 0
    for cur in np2:
        print(cur)
        index += 1
        if index>20:
            break
