<?php

/**
 * Created by Mr. Jason A. Mullings.
 * File Name: AI_Vectors.php
 * User: jlmconsulting
 * Date: 01/05/2016
 * Time: 10:07 PM
 */
require_once(dirname(dirname(__FILE__)) . "/AI_Tensors/lib/math.php");
require_once(dirname(dirname(__FILE__)) . "/AI_Tensors/AI_Matrices.php");

class AI_Vector
{
    public $dimensions;
    public $values = [];
    public $magnitude = null;

    public function __construct()
    {
        $values = func_get_args();
        if (count($values) == 1 && is_array($values[0])) {
            $values = $values[0];
        }
        foreach ($values as $val) {
            if (!is_numeric($val)) {
                var_dump($val);
                throw new \InvalidArgumentException('AI_Vectors values must be numeric');
            }
        }
        $this->dimensions = count($values);
        $this->values = $values;
    }

    public function getDimensions()
    {
        return $this->dimensions;
    }

    public function get($i)
    {
        return $this->values[$i];
    }

    public function getX()
    {
        return $this->get(0);
    }

    public function getY()
    {
        return $this->get(1);
    }

    public function getZ()
    {
        return $this->get(2);
    }

    public function getW()
    {
        return $this->get(3);
    }

    public function __get($key)
    {
        switch ($key) {
            case 'x':
                return $this->getX();
            case 'y':
                return $this->getY();
            case 'z':
                return $this->getZ();
            case 'w':
                return $this->getW();
            case 'all':
                return $this->values;
            case 'xy':
                return array($this->getX(), $this->getY());
            case 'yz':
                return array($this->getY(), $this->getZ());
            case 'xyz':
                return array($this->getX(), $this->getY(), $this->getZ());
        }
        throw new \RuntimeException('Property "' . $key . '" does not exist.');
    }

    public function getMagnitude()
    {
        $sum = 0;
        for ($i = 0; $i < $this->dimensions; $i++) {
            $sum += pow($this->values[$i], 2);
        }
        return $this->magnitude = sqrt($sum);

    }

    public function pnorm($p)
    {
        $sum = 0;
        for ($i = 0; $i < $this->dimensions; $i++) {
            $sum = pow(abs($this->values[$i]), $p);
        }
        return pow($sum, 1 / $p);
    }

    public function powVector($p)
    {
        $sum = array();
        for ($i = 0; $i < $this->dimensions; $i++) {
            $sum[$i] = pow($this->values[$i], $p);
        }
        return $sum;
    }

    public function chebyshevNorm()
    {
        $max = 0;
        for ($i = 0; $i < $this->dimensions; $i++) {
            $val = abs($this->values[$i]);
            if ($val > $max) {
                $max = $val;
            }
        }
        return $max;
    }

    public function normalize()
    {
        $values = array();
        $magnitude = $this->getMagnitude();
        for ($i = 0; $i < $this->dimensions; $i++) {
            $values[] = $this->values[$i] / $magnitude;
        }
        return new AI_Vector($values);
    }

    public function angleTo(AI_Vector $vector)
    {
        return acos($this->dot($vector) / ($this->getMagnitude() * $vector->getMagnitude()));
    }

    public function add(AI_Vector $vector)
    {
        if ($vector->dimensions != $this->dimensions) {
            throw new \InvalidArgumentException('Passed vector requires ' . $this->dimensions
                . ' dimensions, but ' . $vector->dimensions . ' was given.');
        }
        $values = [];
        for ($i = 0; $i < $this->dimensions; $i++) {
            $values[] = $this->values[$i] + $vector->values[$i];
        }
        return new AI_Vector($values);
    }

    public function subtract(AI_Vector $vector)
    {
        return $this->add($vector->multiply(-1));
    }

    public function multiply($scalar)
    {
        $values = [];
        for ($i = 0; $i < $this->dimensions; $i++) {
            $values[] = $this->values[$i] * $scalar;
        }
        return new AI_Vector($values);
    }

    public function dot(AI_Vector $vector)
    {
        if ($vector->dimensions != $this->dimensions) {
            throw new \InvalidArgumentException('Passed vector requires ' . $this->dimensions
                . ' dimensions, but ' . $vector->dimensions . ' was given.');
        }
        $product = 0;
        for ($i = 0; $i < $this->dimensions; $i++) {
            $product += $this->values[$i] * $vector->values[$i];
        }
        return $product;
    }

    public function cross(AI_Vector $vector)
    {
        if ($this->dimensions !== 3 || $vector->dimensions !== 3) {
            throw new \RuntimeException('Cross multiplication can only be done between two 3D vectors.');
        }
        return new AI_Vector(
            $this->getY() * $vector->getZ() - $this->getZ() * $vector->getY(),
            $this->getZ() * $vector->getX() - $this->getX() * $vector->getZ(),
            $this->getX() * $vector->getY() - $this->getY() * $vector->getX()
        );
    }

    public function getLength()
    {
        $sum = null;
        return sqrt(sqr($this->getX()) + sqr($this->getY()) + sqr($this->getZ()));
        foreach ($this->values as $val) {
            $sum += sqr($val);
        }
        return sqrt($sum);

    }

    public function getDistance(AI_Vector $vector)
    {

        $sum = null;
        for ($i = 0; $i < $this->dimensions; $i++) {
            $sum += pow($this->values[$i] - $vector->values[$i], 2);
        }
        return sqrt($sum);

    }


    public function getUnit()
    {

        $len = $this->getLength();
        return new AI_Vector($this->getX() / $len, $this->getY() / $len, $this->getZ() / $len);

    }

    /**
     * Intersection of two lines
     * Four Vectors are required i.e. line1 line2
     * @param AI_Vector $A2
     * @param AI_Vector $B1
     * @param AI_Vector $B2
     * @return bool
     */
    public function getIntersect(AI_Vector $A2, AI_Vector $B1, AI_Vector $B2)
    {
        $A1 = new AI_Vector($this->values);
        $cross1 = new \AI_Vector($B2->subtract($B1)->values);
        $cross2 = new \AI_Vector($A1->subtract($B1)->values);
        $cross3 = new \AI_Vector($A2->subtract($A1)->values);
        $cross4 = new \AI_Vector($B2->subtract($B1)->values);

        $dot1 = new \AI_Vector($cross1->cross($cross2)->values);
        $dot2 = new \AI_Vector($cross3->cross($cross4)->values);
        $dot3 = new \AI_Vector($cross3->cross($cross2)->values);

        $nA = $dot1->dot($dot2);
        $d = $dot2->dot($dot2);
        $nB = $dot3->dot($dot2);

        $t = new \AI_Vector($cross3->multiply(@($nA / $d))->values);
        $A0 = $t->add($A1);

        $s = new \AI_Vector($cross4->multiply(@($nB / $d))->values);
        $B0 = $s->add($B1);

        if ($cross3->dot($A0->subtract($B0)) == 0 || is_nan($cross3->dot($A0->subtract($B0))) == true)
            return true;
        else
            return false;
    }

    /**
     * @param AI_Matrices $mat
     * @param AI_Vector $o point of rotation
     * @return AI_Vector
     */
    public function transformVertex(AI_Matrices $mat, AI_Vector $o)
    {

        return new AI_Vector(
            $this->getX() * $mat->getElement(0, 0) + $this->getY() * $mat->getElement(1, 0) + $this->getZ() * $mat->getElement(2, 0) + $o->getX(),
            $this->getX() * $mat->getElement(0, 1) + $this->getY() * $mat->getElement(1, 1) + $this->getZ() * $mat->getElement(2, 1) + $o->getY(),
            $this->getX() * $mat->getElement(0, 2) + $this->getY() * $mat->getElement(1, 2) + $this->getZ() * $mat->getElement(2, 2) + $o->getZ()
        );

    }

    public function __toString()
    {

        return "(" . $this->getX() . "," . $this->getY() . "," . $this->getZ() . ")";

    }


    public function setX($x)
    {

        $this->values[0] = $x;

    }

    public function setY($y)
    {

        $this->values[1] = $y;

    }

    public function setZ($z)
    {

        $this->values[2] = $z;

    }

    public function mult($r)
    {

        $this->values[0] *= $r;
        $this->values[1] *= $r;
        $this->values[2] *= $r;

    }

    public function output()
    {

        echo "(" . $this->getX() . "," . $this->getY() . "," . $this->getZ() . ")\r\n";

    }

    public function arrayXYZ()
    {

        return $this->values;

    }
}
