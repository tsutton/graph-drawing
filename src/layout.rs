#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
}

impl Vector {
    pub fn distance_to(&self, v: &Vector) -> f64 {
        ((self.x - v.x) * (self.x - v.x) + (self.y - v.y) * (self.y - v.y)).sqrt()
    }

    pub fn unit_vector_to(&self, v: &Vector) -> Vector {
        let distance = self.distance_to(v);
        (*v - *self).scaled(1.0 / distance)
    }

    pub fn scaled(&self, factor: f64) -> Vector {
        Vector {
            x: self.x * factor,
            y: self.y * factor,
        }
    }
}

impl std::ops::Add for Vector {
    type Output = Vector;

    fn add(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl std::ops::Sub for Vector {
    type Output = Vector;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl From<(f64, f64)> for Vector {
    fn from(pt: (f64, f64)) -> Self {
        Vector { x: pt.0, y: pt.1 }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn vector_test() {
        struct TestCase {
            first: (f64, f64),
            second: (f64, f64),
            expected_distance: f64,
            expected_unit: (f64, f64),
        }
        fn case(first: (f64, f64), second: (f64, f64), dist: f64, unit: (f64, f64)) -> TestCase {
            TestCase {
                first,
                second,
                expected_distance: dist,
                expected_unit: unit,
            }
        }
        let cases = vec![
            case((0.0, 0.0), (1.0, 0.0), 1.0, (1.0, 0.0)),
            case((0.0, 0.0), (-1.0, 0.0), 1.0, (-1.0, 0.0)),
            case((0.0, 0.0), (2.0, 0.0), 2.0, (1.0, 0.0)),
            // case((0.0, 0.0), (3.0, 4.0), 5.0, (0.6, 0.8)),
        ];
        for (i, case) in cases.iter().enumerate() {
            let left: Vector = case.first.into();
            let right: Vector = case.second.into();
            assert_eq!(left.distance_to(&left), 0.0, "case {i}");
            assert_eq!(right.distance_to(&right), 0.0, "case {i}");
            assert_eq!(left.distance_to(&right), case.expected_distance, "case {i}");
            assert_eq!(right.distance_to(&left), case.expected_distance, "case {i}");
            assert_eq!(
                left.unit_vector_to(&right),
                case.expected_unit.into(),
                "case {i}"
            );
        }
    }
}
