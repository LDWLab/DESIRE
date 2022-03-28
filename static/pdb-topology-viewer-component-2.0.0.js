/**
 * pdb-topology-viewer
 * @version 2.0.0
 * @link https://github.com/PDBeurope/pdb-topology-viewer
 * @license Apache 2.0
 */
/**
 * Copyright 2020-2021 Mandar Deshpande <mandar@ebi.ac.uk>
 * European Bioinformatics Institute (EBI, http://www.ebi.ac.uk/)
 * European Molecular Biology Laboratory (EMBL, http://www.embl.de/)
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, 
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and 
 * limitations under the License.
 */
! function () {
    "use strict";
    "SVGPathSeg" in window || (window.SVGPathSeg = function (t, e, i) {
        this.pathSegType = t, this.pathSegTypeAsLetter = e, this._owningPathSegList = i
    }, window.SVGPathSeg.prototype.classname = "SVGPathSeg", window.SVGPathSeg.PATHSEG_UNKNOWN = 0, window.SVGPathSeg.PATHSEG_CLOSEPATH = 1, window.SVGPathSeg.PATHSEG_MOVETO_ABS = 2, window.SVGPathSeg.PATHSEG_MOVETO_REL = 3, window.SVGPathSeg.PATHSEG_LINETO_ABS = 4, window.SVGPathSeg.PATHSEG_LINETO_REL = 5, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS = 6, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL = 7, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS = 8, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL = 9, window.SVGPathSeg.PATHSEG_ARC_ABS = 10, window.SVGPathSeg.PATHSEG_ARC_REL = 11, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS = 12, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL = 13, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS = 14, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL = 15, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS = 16, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL = 17, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS = 18, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL = 19, window.SVGPathSeg.prototype._segmentChanged = function () {
        this._owningPathSegList && this._owningPathSegList.segmentChanged(this)
    }, window.SVGPathSegClosePath = function (t) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CLOSEPATH, "z", t)
    }, window.SVGPathSegClosePath.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegClosePath.prototype.toString = function () {
        return "[object SVGPathSegClosePath]"
    }, window.SVGPathSegClosePath.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter
    }, window.SVGPathSegClosePath.prototype.clone = function () {
        return new window.SVGPathSegClosePath(void 0)
    }, window.SVGPathSegMovetoAbs = function (t, e, i) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_MOVETO_ABS, "M", t), this._x = e, this._y = i
    }, window.SVGPathSegMovetoAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegMovetoAbs.prototype.toString = function () {
        return "[object SVGPathSegMovetoAbs]"
    }, window.SVGPathSegMovetoAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
    }, window.SVGPathSegMovetoAbs.prototype.clone = function () {
        return new window.SVGPathSegMovetoAbs(void 0, this._x, this._y)
    }, Object.defineProperty(window.SVGPathSegMovetoAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegMovetoAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegMovetoRel = function (t, e, i) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_MOVETO_REL, "m", t), this._x = e, this._y = i
    }, window.SVGPathSegMovetoRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegMovetoRel.prototype.toString = function () {
        return "[object SVGPathSegMovetoRel]"
    }, window.SVGPathSegMovetoRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
    }, window.SVGPathSegMovetoRel.prototype.clone = function () {
        return new window.SVGPathSegMovetoRel(void 0, this._x, this._y)
    }, Object.defineProperty(window.SVGPathSegMovetoRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegMovetoRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegLinetoAbs = function (t, e, i) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_ABS, "L", t), this._x = e, this._y = i
    }, window.SVGPathSegLinetoAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoAbs.prototype.toString = function () {
        return "[object SVGPathSegLinetoAbs]"
    }, window.SVGPathSegLinetoAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
    }, window.SVGPathSegLinetoAbs.prototype.clone = function () {
        return new window.SVGPathSegLinetoAbs(void 0, this._x, this._y)
    }, Object.defineProperty(window.SVGPathSegLinetoAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegLinetoAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegLinetoRel = function (t, e, i) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_REL, "l", t), this._x = e, this._y = i
    }, window.SVGPathSegLinetoRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoRel.prototype.toString = function () {
        return "[object SVGPathSegLinetoRel]"
    }, window.SVGPathSegLinetoRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
    }, window.SVGPathSegLinetoRel.prototype.clone = function () {
        return new window.SVGPathSegLinetoRel(void 0, this._x, this._y)
    }, Object.defineProperty(window.SVGPathSegLinetoRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegLinetoRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoCubicAbs = function (t, e, i, n, r, o, s) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS, "C", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r, this._x2 = o, this._y2 = s
    }, window.SVGPathSegCurvetoCubicAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicAbs.prototype.toString = function () {
        return "[object SVGPathSegCurvetoCubicAbs]"
    }, window.SVGPathSegCurvetoCubicAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoCubicAbs.prototype.clone = function () {
        return new window.SVGPathSegCurvetoCubicAbs(void 0, this._x, this._y, this._x1, this._y1, this._x2, this._y2)
    }, Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "x1", {
        get: function () {
            return this._x1
        },
        set: function (t) {
            this._x1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "y1", {
        get: function () {
            return this._y1
        },
        set: function (t) {
            this._y1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "x2", {
        get: function () {
            return this._x2
        },
        set: function (t) {
            this._x2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "y2", {
        get: function () {
            return this._y2
        },
        set: function (t) {
            this._y2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoCubicRel = function (t, e, i, n, r, o, s) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL, "c", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r, this._x2 = o, this._y2 = s
    }, window.SVGPathSegCurvetoCubicRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicRel.prototype.toString = function () {
        return "[object SVGPathSegCurvetoCubicRel]"
    }, window.SVGPathSegCurvetoCubicRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoCubicRel.prototype.clone = function () {
        return new window.SVGPathSegCurvetoCubicRel(void 0, this._x, this._y, this._x1, this._y1, this._x2, this._y2)
    }, Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "x1", {
        get: function () {
            return this._x1
        },
        set: function (t) {
            this._x1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "y1", {
        get: function () {
            return this._y1
        },
        set: function (t) {
            this._y1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "x2", {
        get: function () {
            return this._x2
        },
        set: function (t) {
            this._x2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "y2", {
        get: function () {
            return this._y2
        },
        set: function (t) {
            this._y2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoQuadraticAbs = function (t, e, i, n, r) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS, "Q", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r
    }, window.SVGPathSegCurvetoQuadraticAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticAbs.prototype.toString = function () {
        return "[object SVGPathSegCurvetoQuadraticAbs]"
    }, window.SVGPathSegCurvetoQuadraticAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoQuadraticAbs.prototype.clone = function () {
        return new window.SVGPathSegCurvetoQuadraticAbs(void 0, this._x, this._y, this._x1, this._y1)
    }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "x1", {
        get: function () {
            return this._x1
        },
        set: function (t) {
            this._x1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "y1", {
        get: function () {
            return this._y1
        },
        set: function (t) {
            this._y1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoQuadraticRel = function (t, e, i, n, r) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL, "q", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r
    }, window.SVGPathSegCurvetoQuadraticRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticRel.prototype.toString = function () {
        return "[object SVGPathSegCurvetoQuadraticRel]"
    }, window.SVGPathSegCurvetoQuadraticRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoQuadraticRel.prototype.clone = function () {
        return new window.SVGPathSegCurvetoQuadraticRel(void 0, this._x, this._y, this._x1, this._y1)
    }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "x1", {
        get: function () {
            return this._x1
        },
        set: function (t) {
            this._x1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "y1", {
        get: function () {
            return this._y1
        },
        set: function (t) {
            this._y1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegArcAbs = function (t, e, i, n, r, o, s, a) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_ARC_ABS, "A", t), this._x = e, this._y = i, this._r1 = n, this._r2 = r, this._angle = o, this._largeArcFlag = s, this._sweepFlag = a
    }, window.SVGPathSegArcAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegArcAbs.prototype.toString = function () {
        return "[object SVGPathSegArcAbs]"
    }, window.SVGPathSegArcAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._r1 + " " + this._r2 + " " + this._angle + " " + (this._largeArcFlag ? "1" : "0") + " " + (this._sweepFlag ? "1" : "0") + " " + this._x + " " + this._y
    }, window.SVGPathSegArcAbs.prototype.clone = function () {
        return new window.SVGPathSegArcAbs(void 0, this._x, this._y, this._r1, this._r2, this._angle, this._largeArcFlag, this._sweepFlag)
    }, Object.defineProperty(window.SVGPathSegArcAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "r1", {
        get: function () {
            return this._r1
        },
        set: function (t) {
            this._r1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "r2", {
        get: function () {
            return this._r2
        },
        set: function (t) {
            this._r2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "angle", {
        get: function () {
            return this._angle
        },
        set: function (t) {
            this._angle = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "largeArcFlag", {
        get: function () {
            return this._largeArcFlag
        },
        set: function (t) {
            this._largeArcFlag = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "sweepFlag", {
        get: function () {
            return this._sweepFlag
        },
        set: function (t) {
            this._sweepFlag = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegArcRel = function (t, e, i, n, r, o, s, a) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_ARC_REL, "a", t), this._x = e, this._y = i, this._r1 = n, this._r2 = r, this._angle = o, this._largeArcFlag = s, this._sweepFlag = a
    }, window.SVGPathSegArcRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegArcRel.prototype.toString = function () {
        return "[object SVGPathSegArcRel]"
    }, window.SVGPathSegArcRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._r1 + " " + this._r2 + " " + this._angle + " " + (this._largeArcFlag ? "1" : "0") + " " + (this._sweepFlag ? "1" : "0") + " " + this._x + " " + this._y
    }, window.SVGPathSegArcRel.prototype.clone = function () {
        return new window.SVGPathSegArcRel(void 0, this._x, this._y, this._r1, this._r2, this._angle, this._largeArcFlag, this._sweepFlag)
    }, Object.defineProperty(window.SVGPathSegArcRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "r1", {
        get: function () {
            return this._r1
        },
        set: function (t) {
            this._r1 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "r2", {
        get: function () {
            return this._r2
        },
        set: function (t) {
            this._r2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "angle", {
        get: function () {
            return this._angle
        },
        set: function (t) {
            this._angle = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "largeArcFlag", {
        get: function () {
            return this._largeArcFlag
        },
        set: function (t) {
            this._largeArcFlag = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "sweepFlag", {
        get: function () {
            return this._sweepFlag
        },
        set: function (t) {
            this._sweepFlag = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegLinetoHorizontalAbs = function (t, e) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS, "H", t), this._x = e
    }, window.SVGPathSegLinetoHorizontalAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoHorizontalAbs.prototype.toString = function () {
        return "[object SVGPathSegLinetoHorizontalAbs]"
    }, window.SVGPathSegLinetoHorizontalAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x
    }, window.SVGPathSegLinetoHorizontalAbs.prototype.clone = function () {
        return new window.SVGPathSegLinetoHorizontalAbs(void 0, this._x)
    }, Object.defineProperty(window.SVGPathSegLinetoHorizontalAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegLinetoHorizontalRel = function (t, e) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL, "h", t), this._x = e
    }, window.SVGPathSegLinetoHorizontalRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoHorizontalRel.prototype.toString = function () {
        return "[object SVGPathSegLinetoHorizontalRel]"
    }, window.SVGPathSegLinetoHorizontalRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x
    }, window.SVGPathSegLinetoHorizontalRel.prototype.clone = function () {
        return new window.SVGPathSegLinetoHorizontalRel(void 0, this._x)
    }, Object.defineProperty(window.SVGPathSegLinetoHorizontalRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegLinetoVerticalAbs = function (t, e) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS, "V", t), this._y = e
    }, window.SVGPathSegLinetoVerticalAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoVerticalAbs.prototype.toString = function () {
        return "[object SVGPathSegLinetoVerticalAbs]"
    }, window.SVGPathSegLinetoVerticalAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._y
    }, window.SVGPathSegLinetoVerticalAbs.prototype.clone = function () {
        return new window.SVGPathSegLinetoVerticalAbs(void 0, this._y)
    }, Object.defineProperty(window.SVGPathSegLinetoVerticalAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegLinetoVerticalRel = function (t, e) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL, "v", t), this._y = e
    }, window.SVGPathSegLinetoVerticalRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoVerticalRel.prototype.toString = function () {
        return "[object SVGPathSegLinetoVerticalRel]"
    }, window.SVGPathSegLinetoVerticalRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._y
    }, window.SVGPathSegLinetoVerticalRel.prototype.clone = function () {
        return new window.SVGPathSegLinetoVerticalRel(void 0, this._y)
    }, Object.defineProperty(window.SVGPathSegLinetoVerticalRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoCubicSmoothAbs = function (t, e, i, n, r) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS, "S", t), this._x = e, this._y = i, this._x2 = n, this._y2 = r
    }, window.SVGPathSegCurvetoCubicSmoothAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicSmoothAbs.prototype.toString = function () {
        return "[object SVGPathSegCurvetoCubicSmoothAbs]"
    }, window.SVGPathSegCurvetoCubicSmoothAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoCubicSmoothAbs.prototype.clone = function () {
        return new window.SVGPathSegCurvetoCubicSmoothAbs(void 0, this._x, this._y, this._x2, this._y2)
    }, Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "x2", {
        get: function () {
            return this._x2
        },
        set: function (t) {
            this._x2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "y2", {
        get: function () {
            return this._y2
        },
        set: function (t) {
            this._y2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoCubicSmoothRel = function (t, e, i, n, r) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL, "s", t), this._x = e, this._y = i, this._x2 = n, this._y2 = r
    }, window.SVGPathSegCurvetoCubicSmoothRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicSmoothRel.prototype.toString = function () {
        return "[object SVGPathSegCurvetoCubicSmoothRel]"
    }, window.SVGPathSegCurvetoCubicSmoothRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoCubicSmoothRel.prototype.clone = function () {
        return new window.SVGPathSegCurvetoCubicSmoothRel(void 0, this._x, this._y, this._x2, this._y2)
    }, Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "x2", {
        get: function () {
            return this._x2
        },
        set: function (t) {
            this._x2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "y2", {
        get: function () {
            return this._y2
        },
        set: function (t) {
            this._y2 = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoQuadraticSmoothAbs = function (t, e, i) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS, "T", t), this._x = e, this._y = i
    }, window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype.toString = function () {
        return "[object SVGPathSegCurvetoQuadraticSmoothAbs]"
    }, window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype.clone = function () {
        return new window.SVGPathSegCurvetoQuadraticSmoothAbs(void 0, this._x, this._y)
    }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathSegCurvetoQuadraticSmoothRel = function (t, e, i) {
        window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL, "t", t), this._x = e, this._y = i
    }, window.SVGPathSegCurvetoQuadraticSmoothRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticSmoothRel.prototype.toString = function () {
        return "[object SVGPathSegCurvetoQuadraticSmoothRel]"
    }, window.SVGPathSegCurvetoQuadraticSmoothRel.prototype._asPathString = function () {
        return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
    }, window.SVGPathSegCurvetoQuadraticSmoothRel.prototype.clone = function () {
        return new window.SVGPathSegCurvetoQuadraticSmoothRel(void 0, this._x, this._y)
    }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothRel.prototype, "x", {
        get: function () {
            return this._x
        },
        set: function (t) {
            this._x = t, this._segmentChanged()
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothRel.prototype, "y", {
        get: function () {
            return this._y
        },
        set: function (t) {
            this._y = t, this._segmentChanged()
        },
        enumerable: !0
    }), window.SVGPathElement.prototype.createSVGPathSegClosePath = function () {
        return new window.SVGPathSegClosePath(void 0)
    }, window.SVGPathElement.prototype.createSVGPathSegMovetoAbs = function (t, e) {
        return new window.SVGPathSegMovetoAbs(void 0, t, e)
    }, window.SVGPathElement.prototype.createSVGPathSegMovetoRel = function (t, e) {
        return new window.SVGPathSegMovetoRel(void 0, t, e)
    }, window.SVGPathElement.prototype.createSVGPathSegLinetoAbs = function (t, e) {
        return new window.SVGPathSegLinetoAbs(void 0, t, e)
    }, window.SVGPathElement.prototype.createSVGPathSegLinetoRel = function (t, e) {
        return new window.SVGPathSegLinetoRel(void 0, t, e)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicAbs = function (t, e, i, n, r, o) {
        return new window.SVGPathSegCurvetoCubicAbs(void 0, t, e, i, n, r, o)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicRel = function (t, e, i, n, r, o) {
        return new window.SVGPathSegCurvetoCubicRel(void 0, t, e, i, n, r, o)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticAbs = function (t, e, i, n) {
        return new window.SVGPathSegCurvetoQuadraticAbs(void 0, t, e, i, n)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticRel = function (t, e, i, n) {
        return new window.SVGPathSegCurvetoQuadraticRel(void 0, t, e, i, n)
    }, window.SVGPathElement.prototype.createSVGPathSegArcAbs = function (t, e, i, n, r, o, s) {
        return new window.SVGPathSegArcAbs(void 0, t, e, i, n, r, o, s)
    }, window.SVGPathElement.prototype.createSVGPathSegArcRel = function (t, e, i, n, r, o, s) {
        return new window.SVGPathSegArcRel(void 0, t, e, i, n, r, o, s)
    }, window.SVGPathElement.prototype.createSVGPathSegLinetoHorizontalAbs = function (t) {
        return new window.SVGPathSegLinetoHorizontalAbs(void 0, t)
    }, window.SVGPathElement.prototype.createSVGPathSegLinetoHorizontalRel = function (t) {
        return new window.SVGPathSegLinetoHorizontalRel(void 0, t)
    }, window.SVGPathElement.prototype.createSVGPathSegLinetoVerticalAbs = function (t) {
        return new window.SVGPathSegLinetoVerticalAbs(void 0, t)
    }, window.SVGPathElement.prototype.createSVGPathSegLinetoVerticalRel = function (t) {
        return new window.SVGPathSegLinetoVerticalRel(void 0, t)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicSmoothAbs = function (t, e, i, n) {
        return new window.SVGPathSegCurvetoCubicSmoothAbs(void 0, t, e, i, n)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicSmoothRel = function (t, e, i, n) {
        return new window.SVGPathSegCurvetoCubicSmoothRel(void 0, t, e, i, n)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticSmoothAbs = function (t, e) {
        return new window.SVGPathSegCurvetoQuadraticSmoothAbs(void 0, t, e)
    }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticSmoothRel = function (t, e) {
        return new window.SVGPathSegCurvetoQuadraticSmoothRel(void 0, t, e)
    }, "getPathSegAtLength" in window.SVGPathElement.prototype || (window.SVGPathElement.prototype.getPathSegAtLength = function (t) {
        if (void 0 === t || !isFinite(t)) throw "Invalid arguments.";
        var e = document.createElementNS("http://www.w3.org/2000/svg", "path");
        e.setAttribute("d", this.getAttribute("d"));
        var i = e.pathSegList.numberOfItems - 1;
        if (i <= 0) return 0;
        do {
            if (e.pathSegList.removeItem(i), t > e.getTotalLength()) break;
            i--
        } while (0 < i);
        return i
    })), "SVGPathSegList" in window && "appendItem" in window.SVGPathSegList.prototype || (window.SVGPathSegList = function (t) {
        this._pathElement = t, this._list = this._parsePath(this._pathElement.getAttribute("d")), this._mutationObserverConfig = {
            attributes: !0,
            attributeFilter: ["d"]
        }, this._pathElementMutationObserver = new MutationObserver(this._updateListFromPathMutations.bind(this)), this._pathElementMutationObserver.observe(this._pathElement, this._mutationObserverConfig)
    }, window.SVGPathSegList.prototype.classname = "SVGPathSegList", Object.defineProperty(window.SVGPathSegList.prototype, "numberOfItems", {
        get: function () {
            return this._checkPathSynchronizedToList(), this._list.length
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathElement.prototype, "pathSegList", {
        get: function () {
            return this._pathSegList || (this._pathSegList = new window.SVGPathSegList(this)), this._pathSegList
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathElement.prototype, "normalizedPathSegList", {
        get: function () {
            return this.pathSegList
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathElement.prototype, "animatedPathSegList", {
        get: function () {
            return this.pathSegList
        },
        enumerable: !0
    }), Object.defineProperty(window.SVGPathElement.prototype, "animatedNormalizedPathSegList", {
        get: function () {
            return this.pathSegList
        },
        enumerable: !0
    }), window.SVGPathSegList.prototype._checkPathSynchronizedToList = function () {
        this._updateListFromPathMutations(this._pathElementMutationObserver.takeRecords())
    }, window.SVGPathSegList.prototype._updateListFromPathMutations = function (t) {
        if (this._pathElement) {
            var e = !1;
            t.forEach(function (t) {
                "d" == t.attributeName && (e = !0)
            }), e && (this._list = this._parsePath(this._pathElement.getAttribute("d")))
        }
    }, window.SVGPathSegList.prototype._writeListToPath = function () {
        this._pathElementMutationObserver.disconnect(), this._pathElement.setAttribute("d", window.SVGPathSegList._pathSegArrayAsString(this._list)), this._pathElementMutationObserver.observe(this._pathElement, this._mutationObserverConfig)
    }, window.SVGPathSegList.prototype.segmentChanged = function (t) {
        this._writeListToPath()
    }, window.SVGPathSegList.prototype.clear = function () {
        this._checkPathSynchronizedToList(), this._list.forEach(function (t) {
            t._owningPathSegList = null
        }), this._list = [], this._writeListToPath()
    }, window.SVGPathSegList.prototype.initialize = function (t) {
        return this._checkPathSynchronizedToList(), this._list = [t], (t._owningPathSegList = this)._writeListToPath(), t
    }, window.SVGPathSegList.prototype._checkValidIndex = function (t) {
        if (isNaN(t) || t < 0 || t >= this.numberOfItems) throw "INDEX_SIZE_ERR"
    }, window.SVGPathSegList.prototype.getItem = function (t) {
        return this._checkPathSynchronizedToList(), this._checkValidIndex(t), this._list[t]
    }, window.SVGPathSegList.prototype.insertItemBefore = function (t, e) {
        return this._checkPathSynchronizedToList(), e > this.numberOfItems && (e = this.numberOfItems), t._owningPathSegList && (t = t.clone()), this._list.splice(e, 0, t), (t._owningPathSegList = this)._writeListToPath(), t
    }, window.SVGPathSegList.prototype.replaceItem = function (t, e) {
        return this._checkPathSynchronizedToList(), t._owningPathSegList && (t = t.clone()), this._checkValidIndex(e), ((this._list[e] = t)._owningPathSegList = this)._writeListToPath(), t
    }, window.SVGPathSegList.prototype.removeItem = function (t) {
        this._checkPathSynchronizedToList(), this._checkValidIndex(t);
        var e = this._list[t];
        return this._list.splice(t, 1), this._writeListToPath(), e
    }, window.SVGPathSegList.prototype.appendItem = function (t) {
        return this._checkPathSynchronizedToList(), t._owningPathSegList && (t = t.clone()), this._list.push(t), (t._owningPathSegList = this)._writeListToPath(), t
    }, window.SVGPathSegList._pathSegArrayAsString = function (t) {
        var e = "",
            i = !0;
        return t.forEach(function (t) {
            i ? (i = !1, e += t._asPathString()) : e += " " + t._asPathString()
        }), e
    }, window.SVGPathSegList.prototype._parsePath = function (t) {
        if (!t || 0 == t.length) return [];

        function e() {
            this.pathSegList = []
        }
        var n = this;

        function i(t) {
            this._string = t, this._currentIndex = 0, this._endIndex = this._string.length, this._previousCommand = window.SVGPathSeg.PATHSEG_UNKNOWN, this._skipOptionalSpaces()
        }
        e.prototype.appendSegment = function (t) {
            this.pathSegList.push(t)
        }, i.prototype._isCurrentSpace = function () {
            var t = this._string[this._currentIndex];
            return t <= " " && (" " == t || "\n" == t || "\t" == t || "\r" == t || "\f" == t)
        }, i.prototype._skipOptionalSpaces = function () {
            for (; this._currentIndex < this._endIndex && this._isCurrentSpace();) this._currentIndex++;
            return this._currentIndex < this._endIndex
        }, i.prototype._skipOptionalSpacesOrDelimiter = function () {
            return !(this._currentIndex < this._endIndex && !this._isCurrentSpace() && "," != this._string.charAt(this._currentIndex)) && (this._skipOptionalSpaces() && this._currentIndex < this._endIndex && "," == this._string.charAt(this._currentIndex) && (this._currentIndex++, this._skipOptionalSpaces()), this._currentIndex < this._endIndex)
        }, i.prototype.hasMoreData = function () {
            return this._currentIndex < this._endIndex
        }, i.prototype.peekSegmentType = function () {
            var t = this._string[this._currentIndex];
            return this._pathSegTypeFromChar(t)
        }, i.prototype._pathSegTypeFromChar = function (t) {
            switch (t) {
                case "Z":
                case "z":
                    return window.SVGPathSeg.PATHSEG_CLOSEPATH;
                case "M":
                    return window.SVGPathSeg.PATHSEG_MOVETO_ABS;
                case "m":
                    return window.SVGPathSeg.PATHSEG_MOVETO_REL;
                case "L":
                    return window.SVGPathSeg.PATHSEG_LINETO_ABS;
                case "l":
                    return window.SVGPathSeg.PATHSEG_LINETO_REL;
                case "C":
                    return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS;
                case "c":
                    return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL;
                case "Q":
                    return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS;
                case "q":
                    return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL;
                case "A":
                    return window.SVGPathSeg.PATHSEG_ARC_ABS;
                case "a":
                    return window.SVGPathSeg.PATHSEG_ARC_REL;
                case "H":
                    return window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS;
                case "h":
                    return window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL;
                case "V":
                    return window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS;
                case "v":
                    return window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL;
                case "S":
                    return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS;
                case "s":
                    return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL;
                case "T":
                    return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS;
                case "t":
                    return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL;
                default:
                    return window.SVGPathSeg.PATHSEG_UNKNOWN
            }
        }, i.prototype._nextCommandHelper = function (t, e) {
            return ("+" == t || "-" == t || "." == t || "0" <= t && t <= "9") && e != window.SVGPathSeg.PATHSEG_CLOSEPATH ? e == window.SVGPathSeg.PATHSEG_MOVETO_ABS ? window.SVGPathSeg.PATHSEG_LINETO_ABS : e == window.SVGPathSeg.PATHSEG_MOVETO_REL ? window.SVGPathSeg.PATHSEG_LINETO_REL : e : window.SVGPathSeg.PATHSEG_UNKNOWN
        }, i.prototype.initialCommandIsMoveTo = function () {
            if (!this.hasMoreData()) return !0;
            var t = this.peekSegmentType();
            return t == window.SVGPathSeg.PATHSEG_MOVETO_ABS || t == window.SVGPathSeg.PATHSEG_MOVETO_REL
        }, i.prototype._parseNumber = function () {
            var t = 0,
                e = 0,
                i = 1,
                n = 0,
                r = 1,
                o = 1,
                s = this._currentIndex;
            if (this._skipOptionalSpaces(), this._currentIndex < this._endIndex && "+" == this._string.charAt(this._currentIndex) ? this._currentIndex++ : this._currentIndex < this._endIndex && "-" == this._string.charAt(this._currentIndex) && (this._currentIndex++, r = -1), !(this._currentIndex == this._endIndex || (this._string.charAt(this._currentIndex) < "0" || "9" < this._string.charAt(this._currentIndex)) && "." != this._string.charAt(this._currentIndex))) {
                for (var a = this._currentIndex; this._currentIndex < this._endIndex && "0" <= this._string.charAt(this._currentIndex) && this._string.charAt(this._currentIndex) <= "9";) this._currentIndex++;
                if (this._currentIndex != a)
                    for (var h = this._currentIndex - 1, u = 1; a <= h;) e += u * (this._string.charAt(h--) - "0"), u *= 10;
                if (this._currentIndex < this._endIndex && "." == this._string.charAt(this._currentIndex)) {
                    if (this._currentIndex++, this._currentIndex >= this._endIndex || this._string.charAt(this._currentIndex) < "0" || "9" < this._string.charAt(this._currentIndex)) return;
                    for (; this._currentIndex < this._endIndex && "0" <= this._string.charAt(this._currentIndex) && this._string.charAt(this._currentIndex) <= "9";) i *= 10, n += (this._string.charAt(this._currentIndex) - "0") / i, this._currentIndex += 1
                }
                if (this._currentIndex != s && this._currentIndex + 1 < this._endIndex && ("e" == this._string.charAt(this._currentIndex) || "E" == this._string.charAt(this._currentIndex)) && "x" != this._string.charAt(this._currentIndex + 1) && "m" != this._string.charAt(this._currentIndex + 1)) {
                    if (this._currentIndex++, "+" == this._string.charAt(this._currentIndex) ? this._currentIndex++ : "-" == this._string.charAt(this._currentIndex) && (this._currentIndex++, o = -1), this._currentIndex >= this._endIndex || this._string.charAt(this._currentIndex) < "0" || "9" < this._string.charAt(this._currentIndex)) return;
                    for (; this._currentIndex < this._endIndex && "0" <= this._string.charAt(this._currentIndex) && this._string.charAt(this._currentIndex) <= "9";) t *= 10, t += this._string.charAt(this._currentIndex) - "0", this._currentIndex++
                }
                var c = e + n;
                if (c *= r, t && (c *= Math.pow(10, o * t)), s != this._currentIndex) return this._skipOptionalSpacesOrDelimiter(), c
            }
        }, i.prototype._parseArcFlag = function () {
            if (!(this._currentIndex >= this._endIndex)) {
                var t = !1,
                    e = this._string.charAt(this._currentIndex++);
                if ("0" == e) t = !1;
                else {
                    if ("1" != e) return;
                    t = !0
                }
                return this._skipOptionalSpacesOrDelimiter(), t
            }
        }, i.prototype.parseSegment = function () {
            var t = this._string[this._currentIndex],
                e = this._pathSegTypeFromChar(t);
            if (e == window.SVGPathSeg.PATHSEG_UNKNOWN) {
                if (this._previousCommand == window.SVGPathSeg.PATHSEG_UNKNOWN) return null;
                if ((e = this._nextCommandHelper(t, this._previousCommand)) == window.SVGPathSeg.PATHSEG_UNKNOWN) return null
            } else this._currentIndex++;
            switch (this._previousCommand = e) {
                case window.SVGPathSeg.PATHSEG_MOVETO_REL:
                    return new window.SVGPathSegMovetoRel(n, this._parseNumber(), this._parseNumber());
                case window.SVGPathSeg.PATHSEG_MOVETO_ABS:
                    return new window.SVGPathSegMovetoAbs(n, this._parseNumber(), this._parseNumber());
                case window.SVGPathSeg.PATHSEG_LINETO_REL:
                    return new window.SVGPathSegLinetoRel(n, this._parseNumber(), this._parseNumber());
                case window.SVGPathSeg.PATHSEG_LINETO_ABS:
                    return new window.SVGPathSegLinetoAbs(n, this._parseNumber(), this._parseNumber());
                case window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL:
                    return new window.SVGPathSegLinetoHorizontalRel(n, this._parseNumber());
                case window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS:
                    return new window.SVGPathSegLinetoHorizontalAbs(n, this._parseNumber());
                case window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL:
                    return new window.SVGPathSegLinetoVerticalRel(n, this._parseNumber());
                case window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS:
                    return new window.SVGPathSegLinetoVerticalAbs(n, this._parseNumber());
                case window.SVGPathSeg.PATHSEG_CLOSEPATH:
                    return this._skipOptionalSpaces(), new window.SVGPathSegClosePath(n);
                case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL:
                    var i = {
                        x1: this._parseNumber(),
                        y1: this._parseNumber(),
                        x2: this._parseNumber(),
                        y2: this._parseNumber(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    };
                    return new window.SVGPathSegCurvetoCubicRel(n, i.x, i.y, i.x1, i.y1, i.x2, i.y2);
                case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS:
                    return i = {
                        x1: this._parseNumber(),
                        y1: this._parseNumber(),
                        x2: this._parseNumber(),
                        y2: this._parseNumber(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegCurvetoCubicAbs(n, i.x, i.y, i.x1, i.y1, i.x2, i.y2);
                case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL:
                    return i = {
                        x2: this._parseNumber(),
                        y2: this._parseNumber(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegCurvetoCubicSmoothRel(n, i.x, i.y, i.x2, i.y2);
                case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS:
                    return i = {
                        x2: this._parseNumber(),
                        y2: this._parseNumber(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegCurvetoCubicSmoothAbs(n, i.x, i.y, i.x2, i.y2);
                case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL:
                    return i = {
                        x1: this._parseNumber(),
                        y1: this._parseNumber(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegCurvetoQuadraticRel(n, i.x, i.y, i.x1, i.y1);
                case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS:
                    return i = {
                        x1: this._parseNumber(),
                        y1: this._parseNumber(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegCurvetoQuadraticAbs(n, i.x, i.y, i.x1, i.y1);
                case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL:
                    return new window.SVGPathSegCurvetoQuadraticSmoothRel(n, this._parseNumber(), this._parseNumber());
                case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS:
                    return new window.SVGPathSegCurvetoQuadraticSmoothAbs(n, this._parseNumber(), this._parseNumber());
                case window.SVGPathSeg.PATHSEG_ARC_REL:
                    return i = {
                        x1: this._parseNumber(),
                        y1: this._parseNumber(),
                        arcAngle: this._parseNumber(),
                        arcLarge: this._parseArcFlag(),
                        arcSweep: this._parseArcFlag(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegArcRel(n, i.x, i.y, i.x1, i.y1, i.arcAngle, i.arcLarge, i.arcSweep);
                case window.SVGPathSeg.PATHSEG_ARC_ABS:
                    return i = {
                        x1: this._parseNumber(),
                        y1: this._parseNumber(),
                        arcAngle: this._parseNumber(),
                        arcLarge: this._parseArcFlag(),
                        arcSweep: this._parseArcFlag(),
                        x: this._parseNumber(),
                        y: this._parseNumber()
                    }, new window.SVGPathSegArcAbs(n, i.x, i.y, i.x1, i.y1, i.arcAngle, i.arcLarge, i.arcSweep);
                default:
                    throw "Unknown path seg type."
            }
        };
        var r = new e,
            o = new i(t);
        if (!o.initialCommandIsMoveTo()) return [];
        for (; o.hasMoreData();) {
            var s = o.parseSegment();
            if (!s) return [];
            r.appendSegment(s)
        }
        return r.pathSegList
    })
}();
var __awaiter = this && this.__awaiter || function (t, s, a, h) {
        return new(a = a || Promise)(function (i, e) {
            function n(t) {
                try {
                    o(h.next(t))
                } catch (t) {
                    e(t)
                }
            }

            function r(t) {
                try {
                    o(h.throw(t))
                } catch (t) {
                    e(t)
                }
            }

            function o(t) {
                var e;
                t.done ? i(t.value) : ((e = t.value) instanceof a ? e : new a(function (t) {
                    t(e)
                })).then(n, r)
            }
            o((h = h.apply(t, s || [])).next())
        })
    },
    __generator = this && this.__generator || function (i, n) {
        var r, o, s, t, a = {
            label: 0,
            sent: function () {
                if (1 & s[0]) throw s[1];
                return s[1]
            },
            trys: [],
            ops: []
        };
        return t = {
            next: e(0),
            throw: e(1),
            return: e(2)
        }, "function" == typeof Symbol && (t[Symbol.iterator] = function () {
            return this
        }), t;

        function e(e) {
            return function (t) {
                return function (e) {
                    if (r) throw new TypeError("Generator is already executing.");
                    for (; a;) try {
                        if (r = 1, o && (s = 2 & e[0] ? o.return : e[0] ? o.throw || ((s = o.return) && s.call(o), 0) : o.next) && !(s = s.call(o, e[1])).done) return s;
                        switch (o = 0, s && (e = [2 & e[0], s.value]), e[0]) {
                            case 0:
                            case 1:
                                s = e;
                                break;
                            case 4:
                                return a.label++, {
                                    value: e[1],
                                    done: !1
                                };
                            case 5:
                                a.label++, o = e[1], e = [0];
                                continue;
                            case 7:
                                e = a.ops.pop(), a.trys.pop();
                                continue;
                            default:
                                if (!(s = 0 < (s = a.trys).length && s[s.length - 1]) && (6 === e[0] || 2 === e[0])) {
                                    a = 0;
                                    continue
                                }
                                if (3 === e[0] && (!s || e[1] > s[0] && e[1] < s[3])) {
                                    a.label = e[1];
                                    break
                                }
                                if (6 === e[0] && a.label < s[1]) {
                                    a.label = s[1], s = e;
                                    break
                                }
                                if (s && a.label < s[2]) {
                                    a.label = s[2], a.ops.push(e);
                                    break
                                }
                                s[2] && a.ops.pop(), a.trys.pop();
                                continue
                        }
                        e = n.call(i, a)
                    } catch (t) {
                        e = [6, t], o = 0
                    } finally {
                        r = s = 0
                    }
                    if (5 & e[0]) throw e[1];
                    return {
                        value: e[0] ? e[1] : void 0,
                        done: !0
                    }
                }([e, t])
            }
        }
    },
    interpolateLinearly = window.interpolateLinearly,
    RdPu = window.RdPu,
    YlGn = window.YlGn,
    aaPropertyConstants = window.aaPropertyConstants,
    aaColorData = window.aaColorData,
    masking_range_array = window.masking_range_array,
    masked_array = window.masked_array,
    parsePVData = window.parsePVData,
    viewerInstance = window.viewerInstance,
    selectSections_RV1 = window.selectSections_RV1,
    rv3VUEcomponent = window.vm,
    filterRange = window.filterRange ? window.filterRange : "-10000,10000",
    PdbTopologyViewerPlugin = function () {
        function t() {
            this.defaultColours = {
                domainSelection: "rgb(255,0,0)",
                mouseOver: "rgb(105,105,105)",
                borderColor: "rgb(0,0,0)",
                qualityGreen: "rgb(0,182.85714285714286,0)",
                qualityRed: "rgb(291.42857142857144,0,0)",
                qualityYellow: "rgb(364.2857142857143,364.2857142857143,75.71428571428572)",
                qualityRiboVision: "rgb(203,203,203)",
                qualityOrange: "rgb(291.42857142857144,121.42857142857143,0)",
                qualityBlank: "rgb(255,255,255)"
            }, this.displayStyle = "border:1px solid #696969;", this.errorStyle = "border:1px solid #696969; height:54%; padding-top:46%; text-align:center; font-weight:bold;", this.menuStyle = "position:relative;height:38px;line-height:38px;background-color:#96c9dc;padding: 0 10px;font-size:16px; color: black;", this.svgWidth = 100, this.svgHeight = 100, this.pvAPI = !1, this.subscribeEvents = !0, this.createNewEvent = function (t) {
                var n = {};
                return t.forEach(function (t, e) {
                    var i;
                    "function" == typeof MouseEvent ? i = new MouseEvent(t, {
                        view: window,
                        bubbles: !0,
                        cancelable: !0
                    }) : "function" == typeof document.createEvent && (i = document.createEvent("MouseEvents")).initEvent(t, !0, !0), n[t] = i
                }), n
            }, this.getAnnotationFromMappings = function () {
                var r = this;
                if (void 0 !== this.apiData[1]) {
                    function t(t) {
                        if (void 0 !== o[s[t]] && 0 !== Object.entries(o[s[t]]).length) {
                            var e = [],
                                i = o[s[t]];
                            for (var n in i) i[n].mappings.forEach(function (t) {
                                t.entity_id == r.entityId && t.chain_id == r.chainId && e.push({
                                    start: t.start.residue_number,
                                    end: t.end.residue_number,
                                    color: void 0
                                })
                            });
                            0 < e.length && a.domainTypes.push({
                                label: s[t],
                                data: e
                            })
                        }
                    }
                    for (var o = this.apiData[1][this.entryId], s = ["UniProt", "CATH", "Pfam", "SCOP", "RiboVision"], a = this, e = 0; e < 3; e++) t(e)
                }
            }, this.createDomainDropdown = function () {
                if (void 0 === this.domainTypes && (this.domainTypes = [{
                        label: "Select data",
                        data: null
                    }], this.getAnnotationFromRibovision(this.local_mapped_aa_properties), this.selectedDomain = this.domainTypes[0]), 1 < this.domainTypes.length) {
                    var i = "";
                    this.domainTypes.forEach(function (t, e) {
                        i = i + '<option value="' + e + '">' + t.label + "</option>"
                    });
                    var t = this.targetEle.querySelector(".menuSelectbox");
                    t.innerHTML = i, t.addEventListener("change", this.displayDomain.bind(this)), t.addEventListener("change", this.updateProperty.bind(this)), this.targetEle.querySelector(".resetIcon").addEventListener("click", this.resetDisplay.bind(this)), this.targetEle.querySelector(".saveSVG").addEventListener("click", this.saveSVG.bind(this)), rv3VUEcomponent.topology_loaded = !0
                }
            }
        }
        return t.prototype.render = function (t, e) {
            var i = this;
            e && void 0 !== e.displayStyle && null != e.displayStyle && (this.displayStyle += e.displayStyle), e && void 0 !== e.errorStyle && null != e.errorStyle && (this.errorStyle += e.errorStyle), e && void 0 !== e.menuStyle && null != e.menuStyle && (this.menuStyle += e.menuStyle), this.targetEle = t, this.targetEle && (this.targetEle.innerHTML = ""), t && e && e.entryId && e.entityId ? (0 == e.subscribeEvents && (this.subscribeEvents = !1), 1 == e.pvAPI && (this.pvAPI = !0), this.entityId = e.entityId, this.entryId = e.entryId.toLowerCase(), void 0 === e.chainId || null == e.chainId ? this.getObservedResidues(this.entryId).then(function (t) {
                void 0 !== t && void 0 !== t[i.entryId] && void 0 !== t[i.entryId][i.entityId] ? (i.chainId = t[i.entryId][i.entityId][0].chain_id, i.initPainting()) : i.displayError()
            }) : (this.chainId = e.chainId, this.initPainting())) : this.displayError("param")
        }, t.prototype.initPainting = function () {
            var e = this;
            this.alreadyRan || (this.alreadyRan = !0, this.getApiData(this.entryId, this.chainId).then(function (t) {
                if (t) {
                    if (void 0 === t[0] || void 0 === t[2] || void 0 === t[4]) return void e.displayError();
                    e.apiData = t, e.pdbevents = e.createNewEvent(["PDB.topologyViewer" + e.suffix + ".click", "PDB.topologyViewer" + e.suffix + ".mouseover", "PDB.topologyViewer" + e.suffix + ".mouseout"]), e.getPDBSequenceArray(e.apiData[0][e.entryId]), e.drawTopologyStructures(), e.createDomainDropdown(), e.subscribeEvents && e.subscribeWcEvents()
                }
            }))
        }, t.prototype.displayError = function (t) {
            var e = "Error: Data not available!";
            "param" == t && (e = "Error: Invalid Parameters!"), this.targetEle && (this.targetEle.innerHTML = '<div style="' + this.errorStyle + '">' + e + "</div>")
        }, t.prototype.getObservedResidues = function (i) {
            return __awaiter(this, void 0, void 0, function () {
                var e;
                return __generator(this, function (t) {
                    switch (t.label) {
                        case 0:
                            return t.trys.push([0, 3, , 4]), [4, fetch("https://www.ebi.ac.uk/pdbe/api/pdb/entry/observed_residues_ratio/" + i)];
                        case 1:
                            return [4, t.sent().json()];
                        case 2:
                            return [2, t.sent()];
                        case 3:
                            return e = t.sent(), console.log("Couldn't load UniProt variants", e), [3, 4];
                        case 4:
                            return [2]
                    }
                })
            })
        }, t.prototype.getApiData = function (i, n) {
            return __awaiter(this, void 0, void 0, function () {
                var e;
                return __generator(this, function (t) {
                    return e = ["https://www.ebi.ac.uk/pdbe/api/pdb/entry/entities/" + i, "https://www.ebi.ac.uk/pdbe/api/mappings/" + i, "https://www.ebi.ac.uk/pdbe/api/topology/entry/" + i, "https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/" + i, "https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/" + i + "/chain/" + n], this.pvAPI && (e = ["/proOrigamiTopology/ENTITY-" + i + "-" + this.entityId + "-" + n, "/proOrigamiTopology/EMPTY", "/proOrigamiTopology/TOPOLOGY-" + i + "-" + this.entityId + "-" + n, "/proOrigamiTopology/EMPTY", "/proOrigamiTopology/COVERAGE-" + i + "-" + this.entityId + "-" + n]), [2, Promise.all(e.map(function (t) {
                        return fetch(t)
                    })).then(function (t) {
                        return Promise.all(t.map(function (t) {
                            return 200 == t.status ? t.json() : void 0
                        }))
                    })]
                })
            })
        }, t.prototype.getPDBSequenceArray = function (t) {
            for (var e = t.length, i = 0; i < e; i++) t[i].entity_id == this.entityId && (this.sequenceArr = t[i].sequence.split(""))
        }, t.prototype.chunkArray = function (t, e) {
            for (var i = [], n = 0, r = t.length; n < r;) i.push(t.slice(n, n += e));
            return i
        }, t.prototype.getDomainRange = function () {
            var e = this,
                i = [],
                t = this.apiData[2][this.entryId][this.entityId][this.chainId],
                n = Number(filterRange.split(",")[0]),
                r = Number(filterRange.split(",")[1]);
            if (5e3 < r) {
                for (var o = this.apiData[2][this.entryId][this.entityId][this.chainId].coils, s = this.apiData[2][this.entryId][this.entityId][this.chainId].helices, a = this.apiData[2][this.entryId][this.entityId][this.chainId].strands, h = 0; h < s.length; h++)(s[h].start < n && s[h].stop < n || s[h].start > r && s[h].stop > r) && (this.apiData[2][this.entryId][this.entityId][this.chainId].helices.splice(h, 1), h--);
                for (h = 0; h < a.length; h++)(a[h].stop < n || a[h].start > r) && (this.apiData[2][this.entryId][this.entityId][this.chainId].strands.splice(h, 1), h--);
                for (h = 0; h < o.length; h++)(o[h].stop < n && -1 != o[h].start || o[h].start > r && o[h].stop > r) && (this.apiData[2][this.entryId][this.entityId][this.chainId].coils.splice(h, 1), h--);
                for (o = this.apiData[2][this.entryId][this.entityId][this.chainId].coils, h = 0; h < o.length; h++) - 1 == o[0].start && (this.apiData[2][this.entryId][this.entityId][this.chainId].coils.splice(0, 1), h--);
                for (h = o.length - 1; 0 <= h; h--) - 1 == o[o.length - 1].start && (this.apiData[2][this.entryId][this.entityId][this.chainId].coils.splice(o.length - 1, 1), h--);
                for (o = this.apiData[2][this.entryId][this.entityId][this.chainId].coils, h = 0; h < 1; h++)
                    if (o[0].start < n && o[0].stop >= n) {
                        var u = n - o[0].start - 1,
                            c = u;
                        2 < u && (c = 2), this.apiData[2][this.entryId][this.entityId][this.chainId].coils[0].path.splice(0, 2 * (c - 1)), this.apiData[2][this.entryId][this.entityId][this.chainId].coils[0].start = this.apiData[2][this.entryId][this.entityId][this.chainId].coils[0].start + u + 1
                    } for (h = o.length - 1; h > o.length - 2; h--)
                    if (o[o.length - 1].start <= r && o[o.length - 1].stop >= r) {
                        2 < (c = u = o[o.length - 1].stop - 1 - r) && (c = 2);
                        var d = o[o.length - 1].path.length;
                        console.log(u, d), this.apiData[2][this.entryId][this.entityId][this.chainId].coils[o.length - 1].path.splice(d - 2 * (c - 1) - 1, 2 * (c - 1)), this.apiData[2][this.entryId][this.entityId][this.chainId].coils[o.length - 1].stop = this.apiData[2][this.entryId][this.entityId][this.chainId].coils[o.length - 1].stop - u - 1
                    }
            }
            for (var l in t) t[l] && t[l].forEach(function (t) {
                void 0 !== t.path && 0 < t.path.length && (i = i.concat(e.chunkArray(t.path, 2)))
            });
            this.xScale = d3.scaleLinear().domain([d3.min(i, function (t) {
                return t[0]
            }), d3.max(i, function (t) {
                return t[0]
            })]).range([1, this.svgWidth - 1]), this.yScale = d3.scaleLinear().domain([d3.min(i, function (t) {
                return t[1]
            }), d3.max(i, function (t) {
                return t[1]
            })]).range([1, this.svgHeight - 1]), this.zoom = d3.zoom().on("zoom", function () {
                return e.zoomDraw()
            })
        }, t.prototype.drawStrandSubpaths = function (t, e, i) {
            for (var n = this, r = e - t + 1, o = (this.scaledPointsArr[7] - this.scaledPointsArr[1]) / r, s = [], a = 0; a < r; a++) {
                var h = {
                    type: "strands",
                    elementIndex: i
                };
                0 === a ? (h.residue_number = t, h.pathData = [this.scaledPointsArr[4], this.scaledPointsArr[1], this.scaledPointsArr[4], this.scaledPointsArr[1] + o, this.scaledPointsArr[8], this.scaledPointsArr[1] + o, this.scaledPointsArr[8], this.scaledPointsArr[13]]) : (h.residue_number = t + a, h.pathData = [s[a - 1].pathData[2], s[a - 1].pathData[3], s[a - 1].pathData[2], s[a - 1].pathData[3] + o, s[a - 1].pathData[4], s[a - 1].pathData[5] + o, s[a - 1].pathData[4], s[a - 1].pathData[5]]), s.push(h)
            }
            this.svgEle.selectAll(".subpath-strands" + i).remove(), this.svgEle.selectAll(".subpath-strands" + i).data(s).enter().append("path").attr("class", function (t, e) {
                return "strandsSubPath subpath-strands" + i + " topo_res_" + t.residue_number
            }).attr("d", function (t, e) {
                return "M " + t.pathData.join(" ") + " Z"
            }).attr("stroke", "#111").attr("stroke-width", "0").attr("fill", "white").attr("fill-opacity", "0").on("mouseover", function (t) {
                n.mouseoverAction(this, t)
            }).on("mousemove", function (t) {
                n.mouseoverAction(this, t)
            }).on("mouseout", function (t) {
                n.mouseoutAction(this, t)
            }).on("click", function (t) {
                n.clickAction(t)
            })
        }, t.prototype.drawStrandMaskShape = function (i) {
            var n = this.scaledPointsArr,
                t = [],
                e = [];
            this.pvAPI ? (t = [1, 8, 10, 12, 13], e = [0, 2, 3, 4, 5, 7, 9, 11], n[0] > n[6] && (t = [0, 2, 3, 4, 5, 7, 9, 11], e = [1, 8, 10, 12, 13])) : (t = [7, 8, 10, 12], e = [0, 1, 2, 3, 4, 5, 9, 11, 13], n[0] > n[6] && (t = [0, 1, 2, 3, 4, 5, 9, 11, 13], e = [7, 8, 10, 12]));
            for (var r = t.length, o = 0; o < r; o++) n[t[o]] = n[t[o]] + .3;
            for (var s = e.length, a = 0; a < s; a++) n[e[a]] = n[e[a]] - .3;
            n[14] = n[8], n[15] = n[13], n[16] = n[8], n[17] = n[7], n[18] = n[4], n[19] = n[7], n[20] = n[4], n[21] = n[1], this.svgEle.selectAll(".maskpath-strands" + i).remove(), this.svgEle.selectAll(".maskpath-strands" + i).data([n]).enter().append("path").attr("class", function (t, e) {
                return "strandMaskPath maskpath-strands" + i
            }).attr("d", function (t, e) {
                return "M" + n.join(" ") + "Z"
            }).attr("stroke", "#111").attr("stroke-width", .3).attr("fill", "white").attr("stroke-opacity", 0)
        }, t.prototype.renderTooltip = function (t, e) {
            var i = d3.select(".pdbTopologyTooltip");
            if (null == i._groups[0][0] && (i = d3.select("body").append("div").attr("class", "pdbTopologyTooltip").attr("style", "display: none;width: auto;position: absolute;background: #fff;padding: 5px;border: 1px solid #666;border-radius: 5px;box-shadow: 5px 6px 5px 0 rgba(0,0,0,.17);font-size: .9em;color: #555;z-index: 998;")), "show" === e) {
                var n = d3.event.pageX,
                    r = d3.event.pageY,
                    o = "Residue " + t.residue_number + " (" + this.sequenceArr[t.residue_number - 1] + ")";
                void 0 !== t.tooltipMsg && (o = void 0 !== t.tooltipPosition && "postfix" === t.tooltipPosition ? o + " " + t.tooltipMsg : t.tooltipMsg + " " + o), i.html(o).style("display", "block").style("top", r + 15 + "px").style("left", n + 10 + "px")
            } else i.style("display", "none")
        }, t.prototype.dispatchEvent = function (t, e, i) {
            var n = this.targetEle;
            void 0 !== i && (n = i), void 0 !== e && (this.pdbevents[t].eventData = e), n.dispatchEvent(this.pdbevents[t])
        }, t.prototype.clickAction = function (t) {
            1 == masked_array[t.residue_number - 1] && this.dispatchEvent("PDB.topologyViewer" + this.suffix + ".click", {
                residueNumber: t.residue_number,
                type: t.type,
                entryId: this.entryId,
                entityId: this.entityId,
                filterRange: filterRange,
                chainId: this.chainId
            })
        }, t.prototype.mouseoverAction = function (t, e) {
            var i = d3.select(t);
            "NaN" != e.tooltipMsg && (this.renderTooltip(e, "show"), "strands" !== e.type && "helices" !== e.type || i.attr("fill", this.defaultColours.mouseOver).attr("fill-opacity", "0.3"), "coils" === e.type && i.attr("stroke", this.defaultColours.mouseOver).attr("stroke-width", 1), this.dispatchEvent("PDB.topologyViewer" + this.suffix + ".mouseover", {
                residueNumber: e.residue_number,
                type: e.type,
                entryId: this.entryId,
                entityId: this.entityId,
                filterRange: filterRange,
                chainId: this.chainId
            }))
        }, t.prototype.mouseoutAction = function (t, e) {
            var i = "white",
                n = 0,
                r = .3,
                o = d3.select(t);
            this.renderTooltip("", "hide"), o.classed("coloured") ? "coils" === e.type && 0 != masked_array.length && 0 == masked_array[e.residue_number] ? i = this.defaultColours.borderColor : (i = o.attr("data-color"), r = n = 1) : "coils" === e.type && (i = this.defaultColours.borderColor), "strands" !== e.type && "helices" !== e.type || o.attr("fill", i).attr("fill-opacity", n), "coils" === e.type && (o.attr("stroke", i).attr("stroke-opacity", 1), o.attr("stroke", i).attr("stroke-width", r)), this.dispatchEvent("PDB.topologyViewer" + this.suffix + ".mouseout", {
                filterRange: filterRange,
                entryId: this.entryId,
                entityId: this.entityId,
                chainId: this.chainId
            })
        }, t.prototype.drawHelicesSubpaths = function (t, e, i, n) {
            var r = this,
                o = -5;
            this.scaledPointsArr[3] > this.scaledPointsArr[9] && (o = 5);
            var s = e - t + 1;
            o = 0;
            var a = (this.scaledPointsArr[9] - o - this.scaledPointsArr[3]) / s,
                h = 0,
                u = this.svgEle.select(".helices" + i).node().getBBox().height + a / 2,
                c = u / s;
            h = (a = (u - c / 2) / s) - c / 10, this.scaledPointsArr[3] > this.scaledPointsArr[9] && (h = -(u + c));
            for (var d = [], l = {}, p = 0; p < s; p++) l = {
                type: "helices"
            }, 0 === p ? (this.scaledPointsArr[3] > this.scaledPointsArr[9] ? l.residue_number = e : l.residue_number = t, l.pathData = [this.scaledPointsArr[0], this.scaledPointsArr[3] + h, this.scaledPointsArr[4], this.scaledPointsArr[3] + h, this.scaledPointsArr[4], this.scaledPointsArr[3] + h + a, this.scaledPointsArr[0], this.scaledPointsArr[3] + h + a]) : (this.scaledPointsArr[3] > this.scaledPointsArr[9] ? l.residue_number = e - p : l.residue_number = t + p, l.pathData = [d[p - 1].pathData[6], d[p - 1].pathData[7], d[p - 1].pathData[4], d[p - 1].pathData[5], d[p - 1].pathData[4], d[p - 1].pathData[5] + a, d[p - 1].pathData[6], d[p - 1].pathData[5] + a]), d.push(l);
            this.svgEle.selectAll(".subpath-helices" + i).remove(), this.svgEle.selectAll(".subpath-helices" + i).data(d).enter().append("path").attr("class", function (t) {
                return "helicesSubPath subpath-helices" + i + " topo_res_" + t.residue_number
            }).attr("d", function (t) {
                return "M" + t.pathData.join(" ") + " Z"
            }).attr("stroke", "#111").attr("stroke-width", "0").attr("fill", "white").attr("fill-opacity", "0").on("mouseover", function (t) {
                r.mouseoverAction(this, t)
            }).on("mousemove", function (t) {
                r.mouseoverAction(this, t)
            }).on("mouseout", function (t) {
                r.mouseoutAction(this, t)
            }).on("click", function (t) {
                r.clickAction(t)
            })
        }, t.prototype.drawHelicesMaskShape = function (e) {
            var t = [
                [this.scaledPointsArr[0] - .3, this.scaledPointsArr[1], this.scaledPointsArr[2], this.scaledPointsArr[3] - .3, this.scaledPointsArr[4] + .3, this.scaledPointsArr[5], this.scaledPointsArr[4] + .3, this.scaledPointsArr[3], this.scaledPointsArr[0] - .3, this.scaledPointsArr[3]],
                [this.scaledPointsArr[6] + .3, this.scaledPointsArr[7], this.scaledPointsArr[8], this.scaledPointsArr[9] + .3, this.scaledPointsArr[10] - .3, this.scaledPointsArr[11], this.scaledPointsArr[10] - .3, this.scaledPointsArr[9], this.scaledPointsArr[6] + .3, this.scaledPointsArr[9]]
            ];
            this.scaledPointsArr[3] > this.scaledPointsArr[9] && (t = [
                [this.scaledPointsArr[0] - .3, this.scaledPointsArr[1], this.scaledPointsArr[2], this.scaledPointsArr[3] + 2, this.scaledPointsArr[4] + .3, this.scaledPointsArr[5], this.scaledPointsArr[4] + .3, this.scaledPointsArr[3], this.scaledPointsArr[0] - .3, this.scaledPointsArr[3]],
                [this.scaledPointsArr[6] + .3, this.scaledPointsArr[7], this.scaledPointsArr[8], this.scaledPointsArr[9] - .3, this.scaledPointsArr[10] - .3, this.scaledPointsArr[11], this.scaledPointsArr[10] - .3, this.scaledPointsArr[9], this.scaledPointsArr[6] + .3, this.scaledPointsArr[9]]
            ]), this.svgEle.selectAll(".maskpath-helices" + e).remove(), this.svgEle.selectAll(".maskpath-helices" + e).data(t).enter().append("path").attr("class", function (t) {
                return "helicesMaskPath maskpath-helices" + e
            }).attr("d", function (t) {
                return "M" + t[0] + " " + t[1] + " Q" + t[2] + " " + t[3] + " " + t[4] + " " + t[5] + " L" + t[6] + " " + t[7] + " " + t[8] + " " + t[9] + " Z"
            }).attr("stroke", "#111").attr("stroke-width", .3).attr("fill", "white").attr("stroke-opacity", 0)
        }, t.prototype.drawCoilsSubpaths = function (t, e, i) {
            var n = this,
                r = this.svgEle.select(".coils" + i),
                o = e - t + 1,
                s = r.node().getTotalLength() / o,
                a = [],
                h = void 0,
                u = void 0,
                c = {};
            if (1 == o) c = {
                residue_number: t,
                type: "coils",
                pathData: n.scaledPointsArr,
                elementIndex: i
            }, a.push(c);
            else
                for (var d = 0; d < o; d++) {
                    var l = s * (d + 1),
                        p = r.node().getPointAtLength(l),
                        S = r.node().getPathSegAtLength(l);
                    c = {
                        residue_number: t + d,
                        type: "coils",
                        elementIndex: i
                    }, 1 === S ? c.pathData = n.scaledPointsArr.slice(0, 2) : (void 0 === u ? c.pathData = n.scaledPointsArr.slice(0, 2 * S) : (c.pathData = n.scaledPointsArr.slice(2 * u, 2 * S), c.pathData.unshift(h.x, h.y)), h = p, u = S), c.pathData = c.pathData.concat([p.x, p.y]), a.push(c)
                } - 1 !== t && -1 !== e && (this.svgEle.selectAll(".subpath-coils" + i).remove(), this.svgEle.selectAll(".subpath-coils" + i).data(a).enter().append("path").attr("class", function (t) {
                    return "coilsSubPath subpath-coils" + i + " topo_res_" + t.residue_number
                }).attr("d", function (t) {
                    return "M " + t.pathData.join(" ")
                }).attr("stroke", this.defaultColours.borderColor).attr("stroke-width", .3).attr("fill", "none").attr("stroke-opacity", "1").on("mouseover", function (t) {
                    n.mouseoverAction(this, t)
                }).on("mousemove", function (t) {
                    n.mouseoverAction(this, t)
                }).on("mouseout", function (t) {
                    n.mouseoutAction(this, t)
                }).on("click", function (t) {
                    n.clickAction(t)
                }), this.svgEle.selectAll(".coils" + i).attr("stroke-opacity", 0));
            var g = this.apiData[2][this.entryId][this.entityId][this.chainId].terms,
                w = this.apiData[2][this.entryId][this.entityId][this.chainId].coils.length;
            if (0 === i) this.svgEle.selectAll(".terminal_N").remove(), this.svgEle.selectAll(".terminal_N").data([g[0]]).enter().append("text").attr("class", "terminals terminal_N").attr("text-anchor", "middle").text("N").attr("x", a[0].pathData[0]).attr("y", a[0].pathData[1]).attr("stroke", "#0000ff").attr("stroke-width", "0.3").attr("font-size", "3px").attr("style", "-webkit-tap-highlight-color: rgba(0, 0, 0, 0); text-anchor: middle; font-style: normal; font-variant: normal; font-weight: normal; font-stretch: normal; line-height: normal; font-family: Arial;");
            else if (i === w - 1) {
                var _ = a[o - 1].pathData.length,
                    y = -2;
                a[o - 1].pathData[_ - 1] > a[o - 1].pathData[_ - 3] && (y = 2), this.svgEle.selectAll(".terminal_C").remove(), this.svgEle.selectAll(".terminal_C").data([g[1]]).enter().append("text").attr("class", "terminals terminal_C").attr("text-anchor", "middle").text("C").attr("x", a[o - 1].pathData[_ - 2]).attr("y", a[o - 1].pathData[_ - 1] + y).attr("stroke", "#ff0000").attr("stroke-width", "0.3").attr("font-size", "3px").attr("style", "-webkit-tap-highlight-color: rgba(0, 0, 0, 0); text-anchor: middle; font-style: normal; font-variant: normal; font-weight: normal; font-stretch: normal; line-height: normal; font-family: Arial;")
            }
        }, t.prototype.getAdjustedStartAndStop = function (t, e) {
            return null == e ? null : "helices" != t && "coils" != t && "strands" != t || -1 === e.start ? [e.start, e.stop] : e.stop < Number(filterRange.split(",")[0]) || e.start > Number(filterRange.split(",")[1]) ? null : "helices" === t || "strands" === t ? [e.start, e.stop] : e.stop > Number(filterRange.split(",")[1]) ? [e.start, Number(filterRange.split(",")[1])] : e.start < Number(filterRange.split(",")[0]) ? [Number(filterRange.split(",")[0]), e.stop] : [e.start, e.stop]
        }, t.prototype.drawTopologyStructures = function () {
            var c = this;
            this.targetEle.innerHTML = '<div style="' + this.displayStyle + '">\n            <div class="svgSection" style="position:relative;width:100%;"></div>\n            <div style="' + this.menuStyle + '">\n                <a style="color: black;border-bottom:none; cursor:pointer;margin-left: 16px;" target="_blank" href="https://pdbe.org/' + this.entryId + '">' + this.entryId + '</a> | <span class="menuDesc">Entity ' + this.entityId + " | Chain " + this.chainId.toUpperCase() + '</span>\n                <div class="menuOptions" style="float:right;margin-right: 20px;">\n                    <select class="menuSelectbox" style="margin-right: 10px;"><option value="">Select</option></select>\n                    <img class="saveSVG" src="static/alignments/png/Save.png" style="height:15px; width: 15px; border:0;position: relative;margin-right: 15px;cursor:pointer;" title="saveSVG" />\n\n                    <img class="resetIcon" src="static/alignments/png/refresh.png" style="height:15px; width: 15px; border:0;position: absolute;margin-top: 11px;cursor:pointer;" title="Reset view" />\n                </div>\n            </div>\n        </div>';
            var t = this.targetEle.offsetWidth,
                e = this.targetEle.offsetHeight;
            0 == t && (t = this.targetEle.parentElement.offsetWidth), 0 == e && (e = this.targetEle.parentElement.offsetHeight), t <= 330 && (this.targetEle.querySelector(".menuDesc").innerText = this.entityId + " | " + this.chainId.toUpperCase());
            var i = this.targetEle.querySelector(".svgSection"),
                n = e - 40,
                r = t;
            i.style.height = n + "px";
            var o = n - 20,
                s = r - 5;
            i.innerHTML = '<svg class="topoSvg" preserveAspectRatio="xMidYMid meet" viewBox="0 0 100 100" style="width:' + s + "px;height:" + o + 'px;margin:10px 0;"></svg>', this.svgEle = d3.select(this.targetEle).select(".topoSvg"), this.getDomainRange(), this.scaledPointsArr = [], this.svgEle.call(this.zoom).on("contextmenu", function (t, e) {
                d3.event.preventDefault()
            });
            var a = this.apiData[2][this.entryId][this.entityId][this.chainId],
                d = [];

            function h(h) {
                var u = a[h];
                if (!u) return {
                    value: void 0
                };
                u.forEach(function (a, t) {
                    if (void 0 !== a.path && 0 < a.path.length && null != c.getAdjustedStartAndStop(h, a)) {
                        var e = c.getAdjustedStartAndStop(h, a)[0],
                            i = c.getAdjustedStartAndStop(h, a)[1],
                            n = !1;
                        if (-1 == e && -1 == i && "coils" === h) {
                            for (var r = 0, o = 0; o < d.length; o++)(Math.round(a.path[0]) == d[o][0] && (Math.abs(a.path[1] - d[o][1]) < 20 || Math.abs(a.path[1] - d[o][2]) < 20) || Math.round(a.path[a.path.length - 2]) == d[o][0] && (Math.abs(a.path[a.path.length - 1] - d[o][1]) < 20 || Math.abs(a.path[a.path.length - 1] - d[o][2]) < 20)) && (r += 1);
                            2 == r && (n = !0)
                        }
                        if ((-1 !== e || "coils" === h && -1 === e && -1 === i && null != c.getAdjustedStartAndStop(h, u[t + 1]) && null != c.getAdjustedStartAndStop(h, u[t - 1]) || n) && "terms" !== h) {
                            a.secStrType = h, a.pathIndex = t;
                            var s = c.svgEle.selectAll("path." + h + t).data([a]).enter().append("path").attr("class", function () {
                                return -1 === a.start && -1 === a.stop && "terms" !== h ? "dashedEle topologyEle " + h + " " + h + t + " topoEleRange_" + a.start + "-" + a.stop : "topologyEle " + h + " " + h + t + " topoEleRange_" + e + "-" + i
                            }).attr("d", function (t) {
                                for (var e = "M", i = a.path.length, n = !0, r = 0; r < i; r++) {
                                    if ("helices" !== h || 2 !== r && 8 !== r || (e += " Q"), ("helices" === h && 6 === r || "coils" === h && a.path.length < 12 && 8 === r) && (e += " L"), n) {
                                        var o = c.xScale(a.path[r]);
                                        e += " " + o, c.scaledPointsArr.push(o)
                                    } else {
                                        var s = c.yScale(a.path[r]);
                                        e += " " + s, c.scaledPointsArr.push(s)
                                    }
                                    n = !n
                                }
                                return "strands" !== h && "helices" !== h || (e += " Z"), e
                            }).attr("fill", "none").attr("stroke-width", .3).attr("stroke", c.defaultColours.borderColor); - 1 === a.start && -1 === a.stop && s.attr("stroke-dasharray", "0.9"), "strands" === h && (c.drawStrandSubpaths(e, i, t), c.drawStrandMaskShape(t), c.svgEle._groups[0][0].append(s.node())), "helices" === h && (c.drawHelicesSubpaths(e, i, t, 0), c.drawHelicesMaskShape(t), c.svgEle._groups[0][0].append(s.node())), "coils" === h && c.drawCoilsSubpaths(e, i, t), c.scaledPointsArr = []
                        }
                    }
                })
            }
            for (var u in a.helices.forEach(function (t, e) {
                    var i = 0,
                        n = t.path[0] + (t.path[2] - t.path[0]) / 2;
                    i = 1.3 * t.minoraxis * 2, t.path[1] > t.path[3] && (i = 1.3 * t.minoraxis * -2);
                    var r = [t.path[0], t.path[1], n, t.path[1] - i, t.path[2], t.path[1], t.path[2], t.path[3], n, t.path[3] + i, t.path[0], t.path[3]];
                    t.path = r, null != c.getAdjustedStartAndStop("helices", t) && d.push([Math.round(t.path[2]), t.path[3], t.path[6]])
                }), a.strands.forEach(function (t, e) {
                    null != c.getAdjustedStartAndStop("strands", t) && (t.path[9] > t.path[1] ? d.push([Math.round(t.path[6]), t.path[1], t.path[9]]) : d.push([Math.round(t.path[6]), t.path[9], t.path[1]]))
                }), a) {
                var l = h(u);
                if ("object" == typeof l) return l.value
            }
            this.svgEle._groups[0][0].append(this.svgEle.selectAll(".validationResidue").node())
        }, t.prototype.zoomDraw = function () {
            var e = this,
                a = d3.event.transform.rescaleX(this.xScale),
                h = d3.event.transform.rescaleY(this.yScale),
                u = this;
            u.scaledPointsArr = [];
            var t = this.svgEle.selectAll(".topologyEle"),
                c = 0,
                d = 0,
                l = 0;
            t.each(function (t) {
                d3.select(d3.select(this).node()).attr("d", function (t) {
                    c = t.pathIndex, d = t.start, l = t.stop;
                    for (var e = "M", i = t.path.length, n = !0, r = 0; r < i; r++) {
                        if ("helices" !== t.secStrType || 2 !== r && 8 !== r || (e += " Q"), ("helices" === t.secStrType && 6 === r || "coils" === t.secStrType && t.path.length < 12 && 8 === r) && (e += " L"), n) {
                            var o = a(t.path[r]);
                            e += " " + o, u.scaledPointsArr.push(o)
                        } else {
                            var s = h(t.path[r]);
                            e += " " + s, u.scaledPointsArr.push(s)
                        }
                        n = !n
                    }
                    return "strands" !== t.secStrType && "helices" !== t.secStrType || (e += " Z"), e
                }), "helices" === t.secStrType ? (u.drawHelicesSubpaths(d, l, c, 0), u.drawHelicesMaskShape(c), u.svgEle._groups[0][0].append(d3.select(this).node())) : "strands" === t.secStrType ? (u.drawStrandSubpaths(d, l, c), u.drawStrandMaskShape(c), u.svgEle._groups[0][0].append(d3.select(this).node())) : "coils" === t.secStrType && u.drawCoilsSubpaths(d, l, c), u.scaledPointsArr = []
            });
            var s = 0;
            this.svgEle.selectAll(".validationResidue").attr("transform", function (t) {
                var e = u.svgEle.select(".topo_res_" + t.residue_number),
                    i = e.node().getBBox(),
                    n = e.data(),
                    r = {
                        x: 0,
                        y: 0
                    };
                if ("strands" === n[0].type || "helices" === n[0].type) r = {
                    x: i.x + i.width / 2,
                    y: i.y + i.height / 2
                };
                else {
                    var o = e.node().getPointAtLength(e.node().getTotalLength() / 2);
                    r = {
                        x: o.x,
                        y: o.y
                    }
                }
                return s = i.height / 2, "translate(" + r.x + "," + r.y + ")"
            }).attr("d", d3.symbol().type(function (t, e) {
                return d3.symbols[0]
            }).size(s)), this.svgEle.selectAll(".residueSelection").attr("d", function (t) {
                var e = d3.select(this).data();
                return u.svgEle.select(".topo_res_" + e[0].residueNumber).attr("d")
            }), this.svgEle._groups[0][0].querySelectorAll(".coilsSubPath").forEach(function (t) {
                return e.svgEle._groups[0][0].append(t)
            }), this.svgEle._groups[0][0].querySelectorAll(".dashedEle").forEach(function (t) {
                return e.svgEle._groups[0][0].append(t)
            }), this.displayDomain("zoom"), this.svgEle._groups[0][0].querySelectorAll(".validationResidue").forEach(function (t) {
                return e.svgEle._groups[0][0].append(t)
            }), this.svgEle._groups[0][0].querySelectorAll(".residueSelection").forEach(function (t) {
                return e.svgEle._groups[0][0].append(t)
            })
        }, t.prototype.clearHighlight = function () {
            this.svgEle.selectAll(".residueHighlight").remove()
        }, t.prototype.highlight = function (t, e, r, o) {
            function i(e) {
                var t = d.svgEle.select(".topo_res_" + e);
                if (t && t._groups && null == t._groups[0][0]) return {
                    value: void 0
                };
                var i = t.node(),
                    n = t.data();
                r && (a = "string" == typeof r ? h = r : (h = d3.rgb(r.r, r.g, r.b), d3.rgb(r.r, r.g, r.b))), "strands" !== n[0].type && "helices" !== n[0].type ? (a = "none", u = 2, c = .5) : h = "none", d.svgEle.append("path").data([{
                    residueNumber: e
                }]).attr("class", function (t) {
                    return "click" == o ? "residueSelection seletectedResidue_" + e : "residueHighlight highlightResidue_" + e
                }).attr("d", t.attr("d")).attr("fill", a).attr("fill-opacity", .5).attr("stroke", h).attr("stroke-opacity", c).attr("stroke-width", u).on("mouseover", function (t) {
                    s.mouseoverAction(i, n[0])
                }).on("mousemove", function (t) {
                    s.mouseoverAction(i, n[0])
                }).on("mouseout", function (t) {
                    s.mouseoutAction(i, n[0])
                }).on("click", function (t) {
                    s.clickAction(n[0])
                })
            }
            for (var s = this, a = "#000000", h = "#000000", u = .3, c = 0, d = this, n = t; n <= e; n++) {
                var l = i(n);
                if ("object" == typeof l) return l.value
            }
        }, t.prototype.drawValidationShape = function (t, e, i) {
            var n = this,
                r = n.svgEle.select(".topo_res_" + t);
            if (null != r._groups[0][0]) {
                var o = r.node().getBBox(),
                    s = r.data(),
                    a = {
                        x: 0,
                        y: 0
                    };
                if ("strands" === s[0].type || "helices" === s[0].type) a = {
                    x: o.x + o.width / 2,
                    y: o.y + o.height / 2
                };
                else {
                    var h = r.node().getPointAtLength(r.node().getTotalLength() / 2);
                    a = {
                        x: h.x,
                        y: h.y
                    }
                }
                var u = {
                    residue_number: t,
                    tooltipMsg: "Validation issue: RSRZ <br>",
                    tooltipPosition: "prefix"
                };
                this.svgEle.append("path").attr("class", "validationResidue rsrz_" + t).data([u]).attr("fill", i).attr("stroke", "#000").attr("stroke-width", .3).attr("transform", function (t) {
                    return "translate(" + a.x + "," + a.y + ")"
                }).attr("d", d3.symbol().type(function (t, e) {
                    return d3.symbols[0]
                }).size(o.height / 2)).style("display", "none").on("mouseover", function (t) {
                    n.mouseoverAction(this, t)
                }).on("mousemove", function (t) {
                    n.mouseoverAction(this, t)
                }).on("mouseout", function (t) {
                    n.mouseoutAction(this, t)
                }).on("click", function (t) {
                    n.clickAction(t)
                })
            }
        }, t.prototype.getChainStartAndEnd = function () {
            if (this.apiData) {
                for (var t = this.apiData[4][this.entryId].molecules[0].chains, i = {
                        start: 0,
                        end: 0
                    }, e = t.length, n = 0; n < e; n++)
                    if (t[n].chain_id == this.chainId) {
                        t[n].observed.forEach(function (t, e) {
                            0 == e ? (i.start = t.start.residue_number, i.end = t.end.residue_number) : (t.start.residue_number < i.start && (i.start = t.start.residue_number), t.end.residue_number > i.end && (i.end = t.end.residue_number))
                        });
                        break
                    } return i
            }
        }, t.prototype.create2D3DAnnotations = function (r, o, s, t, a, h) {
            var u = this;
            if (t.forEach(function (t, e) {
                    if (a <= e && e <= h) {
                        var i = s.get(e);
                        selectSections_RV1.get(r).push({
                            entity_id: u.entityId,
                            start_residue_number: e,
                            end_residue_number: e,
                            color: i[1],
                            sideChain: !1
                        }), u.defaultColours.qualityRiboVision = "rgb(" + String(i[0].join(",")) + ")";
                        var n = "rgb(" + String(i[0].join(",")) + ")";
                        o.push({
                            start: e,
                            end: e,
                            color: n,
                            tooltipMsg: Number.parseFloat(t).toPrecision(3),
                            tooltipPosition: "prefix"
                        })
                    }
                }), t.size < window.mapped_aa_properties.get("Charge").length)
                for (var e = t.size - 1; e < window.mapped_aa_properties.get("Charge").length; e++) selectSections_RV1.get(r).push({
                    entity_id: u.entityId,
                    start_residue_number: e,
                    end_residue_number: e,
                    color: {
                        r: 255,
                        g: 255,
                        b: 255
                    },
                    sideChain: !1
                });
            return o
        }, t.prototype.getAnnotationFromRibovision = function (t) {
            var l = this,
                p = this.getChainStartAndEnd();
            t && t.forEach(function (t, e) {
                var i = [{}],
                    n = e,
                    r = t;
                selectSections_RV1.set(n, []);
                var aaPropertyConstantsN = aaPropertyConstants.get(n);
                if (aaPropertyConstantsN) {
                    var
                        o = Math.min.apply(Math, aaPropertyConstantsN),
                        s = Math.max.apply(Math, aaPropertyConstantsN),
                        a = aaColorData.get(n),
                        h = parsePVData(r, o, s, a),
                        u = h[0],
                        c = h[1];
                } else {
                    u = new Map();
                    c = new Map();
                    for (let index = p.start; index <= p.end; index++) {
                        let red = 255;
                        let green = 0;
                        let blue = 0;
                        u.set(index, [
                            [
                                red,
                                green,
                                blue
                            ],
                            {
                                r : red,
                                g : green,
                                b : blue
                            }
                        ]);
                        c.set(index, 0);
                    }
                }
                if (selectSections_RV1.get(n).push({
                        entity_id: l.entityId,
                        focus: !0
                    }), void 0 !== c && 0 < (i = l.create2D3DAnnotations(n, i, u, c, p.start, p.end)).length) {
                    var d = l.domainTypes.filter(function (t) {
                        return t.label === n
                    })[0];
                    d && null != d ? d.data = i : l.domainTypes.push({
                        label: n,
                        data: i
                    })
                }
            })
        }, t.prototype.getAnnotationFromOutliers = function () {
            var e = this,
                o = this,
                t = this.getChainStartAndEnd(),
                s = [{
                    start: t.start,
                    end: t.end,
                    color: o.defaultColours.qualityGreen,
                    tooltipMsg: "No validation issue reported for "
                }],
                a = [],
                h = [0];
            if (void 0 !== this.apiData[3]) {
                var i = this.apiData[3][this.entryId];
                void 0 !== i && void 0 !== i.molecules && 0 < i.molecules.length && (i.molecules.forEach(function (t) {
                    t.entity_id == e.entityId && t.chains.forEach(function (t) {
                        t.chain_id == e.chainId && t.models.forEach(function (t) {
                            t.residues.forEach(function (t) {
                                var e = o.defaultColours.qualityYellow,
                                    i = "issue";
                                if (1 !== t.outlier_types.length || "RSRZ" !== t.outlier_types[0]) {
                                    1 === t.outlier_types.length ? e = o.defaultColours.qualityYellow : (e = 2 === t.outlier_types.length ? o.defaultColours.qualityOrange : o.defaultColours.qualityRed, i = "issues"), h.push(t.residue_number);
                                    var n = "Validation " + i + ": " + t.outlier_types.join(", ") + "<br>"; - 1 < a.indexOf(t.residue_number) && (n = "Validation issues: " + t.outlier_types.join(", ") + ", RSRZ<br>"), s.push({
                                        start: parseInt(t.residue_number),
                                        end: parseInt(t.residue_number),
                                        color: e,
                                        tooltipMsg: n,
                                        tooltipPosition: "prefix"
                                    })
                                } else {
                                    e = o.defaultColours.qualityRed, o.drawValidationShape(t.residue_number, "circle", e), a.push(t.residue_number);
                                    var r = h.indexOf(t.residue_number); - 1 < r ? s[r].tooltipMsg = s[r].tooltipMsg.replace("<br>", ", RSRZ<br>") : (s.push({
                                        start: parseInt(t.residue_number),
                                        end: parseInt(t.residue_number),
                                        color: o.defaultColours.qualityGreen,
                                        tooltipMsg: "Validation issue: RSRZ <br>",
                                        tooltipPosition: "prefix"
                                    }), h.push(t.residue_number))
                                }
                            })
                        })
                    })
                }), 0 < s.length && this.domainTypes.push({
                    label: "Quality",
                    data: s
                }))
            }
        }, t.prototype.updateProperty = function () {
            var t = this.targetEle.querySelector(".menuSelectbox"),
                e = parseInt(t.selectedIndex);
            rv3VUEcomponent[this.selected_property_variable_name] = this.domainTypes[e].label
        }, t.prototype.resetTheme = function () {
            var o = this;
            this.svgEle.selectAll(".coloured").each(function (t) {
                var e = d3.select(this),
                    i = e.node();
                e.data()[0].tooltipMsg = void 0, e.data()[0].tooltipPosition = void 0;
                var n = d3.select(i).classed("coloured", !1),
                    r = n.attr("class").split(" "); - 1 < r.indexOf("strandsSubPath") || -1 < r.indexOf("helicesSubPath") ? n.attr("fill", "white").attr("fill-opacity", 0) : n.attr("stroke", o.defaultColours.borderColor).attr("stroke-width", .3)
            }), this.svgEle.selectAll(".validationResidue").style("display", "none")
        }, t.prototype.changeResidueColor = function (e, i, t, n) {
            void 0 === i && (i = this.defaultColours.domainSelection);
            var r = this.svgEle.select(".topo_res_" + e);
            null != r._groups[0][0] && (r.data()[0].tooltipMsg = t, r.data()[0].tooltipPosition = n, r.attr("stroke", function (t) {
                return "coils" === t.type ? i : "#111"
            }).attr("stroke-width", function (t) {
                return "coils" === t.type && !1 === masked_array[e] ? .3 : "coils" === t.type ? 1 : 0
            }).attr("fill", function (t) {
                return "coils" === t.type ? "none" : i
            }).attr("fill-opacity", function (t) {
                return "coils" === t.type ? 0 : 1
            }).classed("coloured", !0).attr("data-color", i))
        }, t.prototype.updateTheme = function (t) {
            var i = this;
            t.forEach(function (t) {
                for (var e = t.start; e <= t.end; e++) i.changeResidueColor(e, t.color, t.tooltipMsg, t.tooltipPosition)
            })
        }, t.prototype.saveSVG = function () {
            var t = this.targetEle.querySelector(".topoSvg"),
                e = this.targetEle.querySelector(".svgSection"),
                i = t.cloneNode(!0),
                n = function (t) {
                    for (var e in void(t = document.createElementNS("http://www.w3.org/2000/svg", t))) t.setAttributeNS(null, e.replace(/[0-9]/g, "o").replace(/\$/g, "d").replace(/\[/g, "b").replace(/[A-Z]/g, function (t, e, i, n) {
                        return "-" + t.toLowerCase()
                    }), (void 0)[e]);
                    return t
                }("svg");
            n.appendChild(i), t && e && e.appendChild(t),
                function (t) {
                    t.setAttribute("xmlns", "http://www.w3.org/2000/svg");
                    var e = t.outerHTML,
                        i = new Blob(['<?xml version="1.0" standalone="no"?>\r\n', e], {
                            type: "image/svg+xml;charset=utf-8"
                        }),
                        n = URL.createObjectURL(i),
                        r = document.createElement("a");
                    r.href = n, r.download = "rv3Topology.svg", document.body.appendChild(r), r.click(), document.body.removeChild(r)
                }(n)
        }, t.prototype.displayDomain = function (t) {
            var e = this.targetEle.querySelector(".menuSelectbox"),
                i = parseInt(e.value);
            if (i) {
                var n = this.domainTypes[i],
                    r = Array.from(aaPropertyConstants.keys());
                if (null !== n.data) {
                    if (this.resetTheme(), this.updateTheme(n.data), r.includes(n.label) && "zoom" !== t) {
                        if ("-10000,10000" != filterRange) var o = filterRange.split(","),
                            s = selectSections_RV1.get(n.label).filter(function (t) {
                                if (t.start_residue_number >= o[0] && t.start_residue_number <= o[1]) return t
                            });
                        else s = selectSections_RV1.get(n.label);
                        viewerInstance.visual.select({
                            data: s,
                            nonSelectedColor: {
                                r: 180,
                                g: 180,
                                b: 180
                            }
                        })
                    }
                    "Quality" === n.label && this.svgEle.selectAll(".validationResidue").style("display", "block")
                } else "zoom" !== t && this.resetTheme()
            } else rv3VUEcomponent.colorSchemeData && this.updateTheme(rv3VUEcomponent.colorSchemeData[0])
        }, t.prototype.resetDisplay = function () {
            this.targetEle.querySelector(".menuSelectbox").value = 0, this.resetTheme(), this.displayDomain()
        }, t.prototype.handleSeqViewerEvents = function (t, e) {
            if (void 0 !== t.eventData) {
                if (t.eventData.entryId.toLowerCase() != this.entryId.toLowerCase() || t.eventData.entityId != this.entityId) return;
                if (t.eventData.elementData.pathData.chain_id && t.eventData.elementData.pathData.chain_id != this.chainId) return;
                var i = "residueSelection";
                "mouseover" == e && (i = "residueHighlight"), this.svgEle.selectAll("." + i).remove();
                var n = void 0,
                    r = void 0;
                if (t.eventData.residueNumber ? (n = t.eventData.residueNumber, r = t.eventData.residueNumber) : t.eventData.elementData.pathData.start.residue_number && t.eventData.elementData.pathData.end.residue_number && (n = t.eventData.elementData.pathData.start.residue_number, r = t.eventData.elementData.pathData.end.residue_number), void 0 !== n && void 0 !== r) {
                    var o;
                    o = t.eventData.elementData.color && 1 == t.eventData.elementData.color.length ? t.eventData.elementData.color[0] : {
                        r: t.eventData.elementData.color[0],
                        g: t.eventData.elementData.color[1],
                        b: t.eventData.elementData.color[2]
                    }, this.highlight(n, r, o, e)
                }
            }
        }, t.prototype.handleProtvistaEvents = function (t, e) {
            if (void 0 !== t.detail) {
                var i = void 0,
                    n = "residueSelection";
                if ("mouseover" == e && (n = "residueHighlight"), this.svgEle.selectAll("." + n).remove(), void 0 !== t.detail.feature) {
                    if (void 0 !== t.detail.feature.accession) {
                        var r = t.detail.feature.accession.split(" ");
                        if ("Chain" == r[0] && r[1].toLowerCase() != this.chainId.toLowerCase()) return
                    } - 1 < t.detail.trackIndex && t.detail.feature.locations && t.detail.feature.locations[0].fragments[t.detail.trackIndex].color && (i = t.detail.feature.locations[0].fragments[t.detail.trackIndex].color), void 0 === i && t.detail.feature.color && (i = t.detail.feature.color)
                }
                void 0 === i && t.detail.color && (i = t.detail.color), void 0 !== i && (i = /rgb/g.test(i) ? i.substring(4, i.length - 1).split(",") : [i]);
                var o = void 0;
                i && (o = 1 == i.length ? t.eventData.elementData.color[0] : {
                    r: t.eventData.elementData.color[0],
                    g: t.eventData.elementData.color[1],
                    b: t.eventData.elementData.color[2]
                }), this.highlight(t.detail.start, t.detail.end, o, e)
            }
        }, t.prototype.handleMolstarEvents = function (t, e) {
            if (void 0 !== t.eventData && 0 < Object.keys(t.eventData).length && (1 == masked_array[t.eventData.seq_id - 1] || null == masked_array[t.eventData.seq_id - 1])) {
                var i = "residueSelection";
                if ("mouseover" == e && (i = "residueHighlight"), this.svgEle.selectAll("." + i).remove(), t.eventData.entry_id.toLowerCase() != this.entryId.toLowerCase() || t.eventData.entity_id != this.entityId) return;
                if (t.eventData.entity_id != this.entityId) return;
                this.highlight(t.eventData.seq_id, t.eventData.seq_id, void 0, e)
            }
        }, t.prototype.subscribeWcEvents = function () {
            var e = this;
            document.addEventListener("PDB.seqViewer.click", function (t) {
                e.handleSeqViewerEvents(t, "click")
            }), document.addEventListener("PDB.seqViewer.mouseover", function (t) {
                e.handleSeqViewerEvents(t, "mouseover")
            }), document.addEventListener("PDB.seqViewer.mouseout", function () {
                e.svgEle.selectAll(".residueHighlight").remove()
            }), document.addEventListener("PDB.litemol.click", function (t) {
                e.svgEle.selectAll(".residueSelection").remove(), t.eventData.entryId.toLowerCase() == e.entryId.toLowerCase() && t.eventData.entityId == e.entityId && t.eventData.chainId.toLowerCase() == e.chainId.toLowerCase() && e.highlight(t.eventData.residueNumber, t.eventData.residueNumber, void 0, "click")
            }), document.addEventListener("PDB.litemol.mouseover", function (t) {
                e.svgEle.selectAll(".residueHighlight").remove(), t.eventData.entryId.toLowerCase() == e.entryId.toLowerCase() && t.eventData.entityId == e.entityId && t.eventData.chainId.toLowerCase() == e.chainId.toLowerCase() && e.highlight(t.eventData.residueNumber, t.eventData.residueNumber, void 0, "mouseover")
            }), document.addEventListener("protvista-click", function (t) {
                e.handleProtvistaEvents(t, "click")
            }), document.addEventListener("protvista-mouseover", function (t) {
                e.handleProtvistaEvents(t, "mouseover")
            }), document.addEventListener("protvista-mouseout", function () {
                e.svgEle.selectAll(".residueHighlight").remove()
            }), document.addEventListener("PDB.molstar.click", function (t) {
                e.handleMolstarEvents(t, "click")
            }), document.addEventListener("PDB.molstar.mouseover", function (t) {
                e.handleMolstarEvents(t, "mouseover")
            }), document.addEventListener("PDB.molstar.mouseout", function () {
                e.svgEle.selectAll(".residueHighlight").remove()
            })
        }, t
    }();
window.PdbTopologyViewerPlugin = PdbTopologyViewerPlugin,
    function (i) {
        var n = {};

        function r(t) {
            if (n[t]) return n[t].exports;
            var e = n[t] = {
                i: t,
                l: !1,
                exports: {}
            };
            return i[t].call(e.exports, e, e.exports, r), e.l = !0, e.exports
        }
        r.m = i, r.c = n, r.d = function (t, e, i) {
            r.o(t, e) || Object.defineProperty(t, e, {
                enumerable: !0,
                get: i
            })
        }, r.r = function (t) {
            "undefined" != typeof Symbol && Symbol.toStringTag && Object.defineProperty(t, Symbol.toStringTag, {
                value: "Module"
            }), Object.defineProperty(t, "__esModule", {
                value: !0
            })
        }, r.t = function (e, t) {
            if (1 & t && (e = r(e)), 8 & t) return e;
            if (4 & t && "object" == typeof e && e && e.__esModule) return e;
            var i = Object.create(null);
            if (r.r(i), Object.defineProperty(i, "default", {
                    enumerable: !0,
                    value: e
                }), 2 & t && "string" != typeof e)
                for (var n in e) r.d(i, n, function (t) {
                    return e[t]
                }.bind(null, n));
            return i
        }, r.n = function (t) {
            var e = t && t.__esModule ? function () {
                return t.default
            } : function () {
                return t
            };
            return r.d(e, "a", e), e
        }, r.o = function (t, e) {
            return Object.prototype.hasOwnProperty.call(t, e)
        }, r.p = "", r(r.s = 7)
    }([function (i, t) {
        function n(t, e) {
            return i.exports = n = Object.setPrototypeOf || function (t, e) {
                return t.__proto__ = e, t
            }, n(t, e)
        }
        i.exports = n
    }, function (e, t) {
        function i(t) {
            return e.exports = i = Object.setPrototypeOf ? Object.getPrototypeOf : function (t) {
                return t.__proto__ || Object.getPrototypeOf(t)
            }, i(t)
        }
        e.exports = i
    }, function (t, e) {
        function n(t, e) {
            for (var i = 0; i < e.length; i++) {
                var n = e[i];
                n.enumerable = n.enumerable || !1, n.configurable = !0, "value" in n && (n.writable = !0), Object.defineProperty(t, n.key, n)
            }
        }
        t.exports = function (t, e, i) {
            return e && n(t.prototype, e), i && n(t, i), t
        }
    }, function (t, e) {
        t.exports = function (t, e) {
            if (!(t instanceof e)) throw new TypeError("Cannot call a class as a function")
        }
    }, function (t, e, i) {
        var n = i(8),
            r = i(9);
        t.exports = function (t, e) {
            return !e || "object" !== n(e) && "function" != typeof e ? r(t) : e
        }
    }, function (t, e, i) {
        var n = i(0);
        t.exports = function (t, e) {
            if ("function" != typeof e && null !== e) throw new TypeError("Super expression must either be null or a function");
            t.prototype = Object.create(e && e.prototype, {
                constructor: {
                    value: t,
                    writable: !0,
                    configurable: !0
                }
            }), e && n(t, e)
        }
    }, function (e, t, i) {
        var n = i(1),
            r = i(0),
            o = i(10),
            s = i(11);

        function a(t) {
            var i = "function" == typeof Map ? new Map : void 0;
            return e.exports = a = function (t) {
                if (null === t || !o(t)) return t;
                if ("function" != typeof t) throw new TypeError("Super expression must either be null or a function");
                if (void 0 !== i) {
                    if (i.has(t)) return i.get(t);
                    i.set(t, e)
                }

                function e() {
                    return s(t, arguments, n(this).constructor)
                }
                return e.prototype = Object.create(t.prototype, {
                    constructor: {
                        value: e,
                        enumerable: !1,
                        writable: !0,
                        configurable: !0
                    }
                }), r(e, t)
            }, a(t)
        }
        e.exports = a
    }, function (t, e, i) {
        "use strict";
        i.r(e);
        var n, r = i(3),
            o = i.n(r),
            s = i(4),
            a = i.n(s),
            h = i(1),
            u = i.n(h),
            c = i(2),
            d = i.n(c),
            l = i(5),
            p = i.n(l),
            S = i(6),
            g = (n = i.n(S)()(HTMLElement), p()(w, n), d()(w, null, [{
                key: "observedAttributes",
                get: function () {
                    return ["pv-aa-properties-variable-name", "pv-selected-property-variable-name", "pv-suffix", "entry-id", "entity-id", "filter-range", "chain-id", "display-style", "error-style", "menu-style", "subscribe-events", "pvapi"]
                }
            }]), d()(w, [{
                key: "validateParams",
                value: function () {
                    return void 0 !== this.entryId && void 0 !== this.entityId && null != this.entryId && null != this.entityId
                }
            }, {
                key: "invokePlugin",
                value: function () {
                    if (this.validateParams()) {
                        void 0 === this.pluginInstance && (this.pluginInstance = new PdbTopologyViewerPlugin);
                        var t = {
                            entryId: this.entryId,
                            entityId: this.entityId,
                            entropyId: this.entropyId
                        };
                        this.pluginInstance.local_mapped_aa_properties = window[this.pv_aa_properties_variable_name];
                        this.pluginInstance.selected_property_variable_name = this.pv_selected_property_variable_name;
                        this.pluginInstance.suffix = this.pv_suffix;
                        void 0 !== this.chainId && null !== this.chainId && (t.chainId = this.chainId), void 0 !== this.displayStyle && null !== this.displayStyle && (t.displayStyle = this.displayStyle), void 0 !== this.errorStyle && null !== this.errorStyle && (t.errorStyle = this.errorStyle), void 0 !== this.menuStyle && null !== this.menuStyle && (t.menuStyle = this.menuStyle), void 0 !== this.subscribeEvents && null !== this.subscribeEvents && (t.subscribeEvents = this.subscribeEvents), void 0 !== this.pvAPI && null !== this.pvAPI && (t.pvAPI = this.pvAPI), void 0 !== this.filterRange && null !== this.filterRange && (t.filterRange = this.filterRange), this.pluginInstance.render(this, t)
                    }
                }
            }, {
                key: "attributeChangedCallback",
                value: function () {
                    this.pv_aa_properties_variable_name = this.getAttribute("pv-aa-properties-variable-name"), this.pv_selected_property_variable_name = this.getAttribute("pv-selected-property-variable-name"), this.pv_suffix = this.getAttribute("pv-suffix"), this.entryId = this.getAttribute("entry-id"), this.entityId = this.getAttribute("entity-id"), this.chainId = this.getAttribute("chain-id"), this.filterRange = this.getAttribute("filter-range"), this.displayStyle = this.getAttribute("display-style"), this.errorStyle = this.getAttribute("error-style"), this.menuStyle = this.getAttribute("menu-style"), this.subscribeEvents = this.getAttribute("subscribe-events"), this.pvAPI = /true/i.test(this.getAttribute("pvapi")), this.invokePlugin()
                }
            }]), w);

        function w() {
            return o()(this, w), a()(this, u()(w).call(this))
        }
        e.default = g, customElements.define("pdb-topology-viewer", g)
    }, function (e, t) {
        function i(t) {
            return "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? e.exports = i = function (t) {
                return typeof t
            } : e.exports = i = function (t) {
                return t && "function" == typeof Symbol && t.constructor === Symbol && t !== Symbol.prototype ? "symbol" : typeof t
            }, i(t)
        }
        e.exports = i
    }, function (t, e) {
        t.exports = function (t) {
            if (void 0 === t) throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
            return t
        }
    }, function (t, e) {
        t.exports = function (t) {
            return -1 !== Function.toString.call(t).indexOf("[native code]")
        }
    }, function (n, t, e) {
        var o = e(0);

        function r(t, e, i) {
            return function () {
                if ("undefined" != typeof Reflect && Reflect.construct && !Reflect.construct.sham) {
                    if ("function" == typeof Proxy) return 1;
                    try {
                        return Date.prototype.toString.call(Reflect.construct(Date, [], function () {})), 1
                    } catch (t) {
                        return
                    }
                }
            }() ? n.exports = r = Reflect.construct : n.exports = r = function (t, e, i) {
                var n = [null];
                n.push.apply(n, e);
                var r = new(Function.bind.apply(t, n));
                return i && o(r, i.prototype), r
            }, r.apply(null, arguments)
        }
        n.exports = r
    }]);