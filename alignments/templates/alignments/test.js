! function(t) {
    function e(r) {
        if (n[r]) return n[r].exports;
        var i = n[r] = {
            exports: {},
            id: r,
            loaded: !1
        };
        return t[r].call(i.exports, i, i.exports, e), i.loaded = !0, i.exports
    }
    var n = {};
    return e.m = t, e.c = n, e.p = "", e(0)
}([function(t, e, n) {
    "use strict";

    function r(t) {
        if (t && t.__esModule) return t;
        var e = {};
        if (null != t)
            for (var n in t) Object.prototype.hasOwnProperty.call(t, n) && (e[n] = t[n]);
        return e["default"] = t, e
    }
    var i = n(111),
        o = r(i);
    n(154);
    var s = o["default"];
    for (var u in o) o.hasOwnProperty(u) && (s[u] = o[u]);
    window && (window.msa = s), t.exports = s
}, function(t, e, n) {
    "use strict";
    t.exports.Model = n(41), t.exports.Collection = n(55), t.exports.Events = n(22), t.exports.extend = n(23)
}, function(t, e, n) {
    "use strict";
    var r = n(60),
        i = n(58),
        o = n(59),
        s = n(5),
        u = function(t) {
            this.cid = r.uniqueId("view"), t || (t = {}), r.extend(this, r.pick(t, l)), this._ensureElement(), this.initialize.apply(this, arguments)
        },
        a = /^(\S+)\s*(.*)$/,
        l = ["model", "collection", "el", "id", "attributes", "className", "tagName", "events"];
    r.extend(u.prototype, i, {
        tagName: "div",
        $: function(t) {
            return this.$el.find(t)
        },
        initialize: function() {},
        render: function() {
            return this
        },
        remove: function() {
            return this._removeElement(), this.stopListening(), this
        },
        _removeElement: function() {
            this.$el.remove()
        },
        setElement: function(t) {
            return this.undelegateEvents(), this._setElement(t), this.delegateEvents(), this
        },
        _setElement: function(t) {
            this.$el = t instanceof s ? t : s(t), this.el = this.$el[0]
        },
        delegateEvents: function(t) {
            if (!t && !(t = r.result(this, "events"))) return this;
            this.undelegateEvents();
            for (var e in t) {
                var n = t[e];
                if (r.isFunction(n) || (n = this[t[e]]), n) {
                    var i = e.match(a);
                    this.delegate(i[1], i[2], r.bind(n, this))
                }
            }
            return this
        },
        delegate: function(t, e, n) {
            this.$el.on(t + ".delegateEvents" + this.cid, e, n)
        },
        undelegateEvents: function() {
            return this.$el && this.$el.off(".delegateEvents" + this.cid), this
        },
        undelegate: function(t, e, n) {
            this.$el.off(t + ".delegateEvents" + this.cid, e, n)
        },
        _createElement: function(t) {
            return document.createElement(t)
        },
        _ensureElement: function() {
            if (this.el) this.setElement(r.result(this, "el"));
            else {
                var t = r.extend({}, r.result(this, "attributes"));
                this.id && (t.id = r.result(this, "id")), this.className && (t["class"] = r.result(this, "className")), this.setElement(this._createElement(r.result(this, "tagName"))), this._setAttributes(t)
            }
        },
        _setAttributes: function(t) {
            this.$el.attr(t)
        }
    }), u.extend = o, t.exports = u
}, function(t, e, n) {
    var r;
    (function(t, i) {
        "use strict";
        var o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        };
        (function() {
            function s(t, e) {
                return t.set(e[0], e[1]), t
            }

            function u(t, e) {
                return t.add(e), t
            }

            function a(t, e, n) {
                switch (n.length) {
                    case 0:
                        return t.call(e);
                    case 1:
                        return t.call(e, n[0]);
                    case 2:
                        return t.call(e, n[0], n[1]);
                    case 3:
                        return t.call(e, n[0], n[1], n[2])
                }
                return t.apply(e, n)
            }

            function l(t, e, n, r) {
                for (var i = -1, o = t ? t.length : 0; ++i < o;) {
                    var s = t[i];
                    e(r, s, n(s), t)
                }
                return r
            }

            function c(t, e) {
                for (var n = -1, r = t ? t.length : 0; ++n < r && e(t[n], n, t) !== !1;);
                return t
            }

            function f(t, e) {
                for (var n = t ? t.length : 0; n-- && e(t[n], n, t) !== !1;);
                return t
            }

            function h(t, e) {
                for (var n = -1, r = t ? t.length : 0; ++n < r;)
                    if (!e(t[n], n, t)) return !1;
                return !0
            }

            function d(t, e) {
                for (var n = -1, r = t ? t.length : 0, i = 0, o = []; ++n < r;) {
                    var s = t[n];
                    e(s, n, t) && (o[i++] = s)
                }
                return o
            }

            function p(t, e) {
                return !!(t ? t.length : 0) && S(t, e, 0) > -1
            }

            function g(t, e, n) {
                for (var r = -1, i = t ? t.length : 0; ++r < i;)
                    if (n(e, t[r])) return !0;
                return !1
            }

            function v(t, e) {
                for (var n = -1, r = t ? t.length : 0, i = Array(r); ++n < r;) i[n] = e(t[n], n, t);
                return i
            }

            function m(t, e) {
                for (var n = -1, r = e.length, i = t.length; ++n < r;) t[i + n] = e[n];
                return t
            }

            function y(t, e, n, r) {
                var i = -1,
                    o = t ? t.length : 0;
                for (r && o && (n = t[++i]); ++i < o;) n = e(n, t[i], i, t);
                return n
            }

            function _(t, e, n, r) {
                var i = t ? t.length : 0;
                for (r && i && (n = t[--i]); i--;) n = e(n, t[i], i, t);
                return n
            }

            function b(t, e) {
                for (var n = -1, r = t ? t.length : 0; ++n < r;)
                    if (e(t[n], n, t)) return !0;
                return !1
            }

            function x(t, e, n) {
                var r;
                return n(t, function(t, n, i) {
                    if (e(t, n, i)) return r = n, !1
                }), r
            }

            function w(t, e, n, r) {
                for (var i = t.length, o = n + (r ? 1 : -1); r ? o-- : ++o < i;)
                    if (e(t[o], o, t)) return o;
                return -1
            }

            function S(t, e, n) {
                if (e !== e) return D(t, n);
                for (var r = n - 1, i = t.length; ++r < i;)
                    if (t[r] === e) return r;
                return -1
            }

            function k(t, e, n, r) {
                for (var i = n - 1, o = t.length; ++i < o;)
                    if (r(t[i], e)) return i;
                return -1
            }

            function j(t, e) {
                var n = t ? t.length : 0;
                return n ? M(t, e) / n : Et
            }

            function O(t, e, n, r, i) {
                return i(t, function(t, i, o) {
                    n = r ? (r = !1, t) : e(n, t, i, o)
                }), n
            }

            function E(t, e) {
                var n = t.length;
                for (t.sort(e); n--;) t[n] = t[n].value;
                return t
            }

            function M(t, e) {
                for (var n, r = -1, i = t.length; ++r < i;) {
                    var o = e(t[r]);
                    o !== Q && (n = n === Q ? o : n + o)
                }
                return n
            }

            function z(t, e) {
                for (var n = -1, r = Array(t); ++n < t;) r[n] = e(n);
                return r
            }

            function A(t, e) {
                return v(e, function(e) {
                    return [e, t[e]]
                })
            }

            function C(t) {
                return function(e) {
                    return t(e)
                }
            }

            function T(t, e) {
                return v(e, function(e) {
                    return t[e]
                })
            }

            function I(t, e) {
                return t.has(e)
            }

            function N(t, e) {
                for (var n = -1, r = t.length; ++n < r && S(e, t[n], 0) > -1;);
                return n
            }

            function L(t, e) {
                for (var n = t.length; n-- && S(e, t[n], 0) > -1;);
                return n
            }

            function R(t) {
                return t && t.Object === Object ? t : null
            }

            function q(t, e) {
                for (var n = t.length, r = 0; n--;) t[n] === e && r++;
                return r
            }

            function F(t) {
                return An[t]
            }

            function P(t) {
                return Cn[t]
            }

            function B(t) {
                return "\\" + In[t]
            }

            function W(t, e) {
                return null == t ? Q : t[e]
            }

            function D(t, e, n) {
                for (var r = t.length, i = e + (n ? 1 : -1); n ? i-- : ++i < r;) {
                    var o = t[i];
                    if (o !== o) return i
                }
                return -1
            }

            function H(t) {
                var e = !1;
                if (null != t && "function" != typeof t.toString) try {
                    e = !!(t + "")
                } catch (n) {}
                return e
            }

            function U(t) {
                for (var e, n = []; !(e = t.next()).done;) n.push(e.value);
                return n
            }

            function V(t) {
                var e = -1,
                    n = Array(t.size);
                return t.forEach(function(t, r) {
                    n[++e] = [r, t]
                }), n
            }

            function $(t, e) {
                for (var n = -1, r = t.length, i = 0, o = []; ++n < r;) {
                    var s = t[n];
                    s !== e && s !== it || (t[n] = it, o[i++] = n)
                }
                return o
            }

            function G(t) {
                var e = -1,
                    n = Array(t.size);
                return t.forEach(function(t) {
                    n[++e] = t
                }), n
            }

            function K(t) {
                var e = -1,
                    n = Array(t.size);
                return t.forEach(function(t) {
                    n[++e] = [t, t]
                }), n
            }

            function X(t) {
                if (!t || !kn.test(t)) return t.length;
                for (var e = wn.lastIndex = 0; wn.test(t);) e++;
                return e
            }

            function Z(t) {
                return t.match(wn)
            }

            function J(t) {
                return Tn[t]
            }

            function Y(t) {
                function e(t) {
                    if (_u(t) && !yf(t) && !(t instanceof i)) {
                        if (t instanceof r) return t;
                        if (Al.call(t, "__wrapped__")) return ho(t)
                    }
                    return new r(t)
                }

                function n() {}

                function r(t, e) {
                    this.__wrapped__ = t, this.__actions__ = [], this.__chain__ = !!e, this.__index__ = 0, this.__values__ = Q
                }

                function i(t) {
                    this.__wrapped__ = t, this.__actions__ = [], this.__dir__ = 1, this.__filtered__ = !1, this.__iteratees__ = [], this.__takeCount__ = Mt, this.__views__ = []
                }

                function R() {
                    var t = new i(this.__wrapped__);
                    return t.__actions__ = ri(this.__actions__), t.__dir__ = this.__dir__, t.__filtered__ = this.__filtered__, t.__iteratees__ = ri(this.__iteratees__), t.__takeCount__ = this.__takeCount__, t.__views__ = ri(this.__views__), t
                }

                function Re() {
                    if (this.__filtered__) {
                        var t = new i(this);
                        t.__dir__ = -1, t.__filtered__ = !0
                    } else t = this.clone(), t.__dir__ *= -1;
                    return t
                }

                function qe() {
                    var t = this.__wrapped__.value(),
                        e = this.__dir__,
                        n = yf(t),
                        r = e < 0,
                        i = n ? t.length : 0,
                        o = Hi(0, i, this.__views__),
                        s = o.start,
                        u = o.end,
                        a = u - s,
                        l = r ? u : s - 1,
                        c = this.__iteratees__,
                        f = c.length,
                        h = 0,
                        d = tc(a, this.__takeCount__);
                    if (!n || i < et || i == a && d == a) return Fr(t, this.__actions__);
                    var p = [];
                    t: for (; a-- && h < d;) {
                        l += e;
                        for (var g = -1, v = t[l]; ++g < f;) {
                            var m = c[g],
                                y = m.iteratee,
                                _ = m.type,
                                b = y(v);
                            if (_ == wt) v = b;
                            else if (!b) {
                                if (_ == xt) continue t;
                                break t
                            }
                        }
                        p[h++] = v
                    }
                    return p
                }

                function Fe(t) {
                    var e = -1,
                        n = t ? t.length : 0;
                    for (this.clear(); ++e < n;) {
                        var r = t[e];
                        this.set(r[0], r[1])
                    }
                }

                function Pe() {
                    this.__data__ = fc ? fc(null) : {}
                }

                function Be(t) {
                    return this.has(t) && delete this.__data__[t]
                }

                function We(t) {
                    var e = this.__data__;
                    if (fc) {
                        var n = e[t];
                        return n === rt ? Q : n
                    }
                    return Al.call(e, t) ? e[t] : Q
                }

                function De(t) {
                    var e = this.__data__;
                    return fc ? e[t] !== Q : Al.call(e, t)
                }

                function He(t, e) {
                    return this.__data__[t] = fc && e === Q ? rt : e, this
                }

                function Ue(t) {
                    var e = -1,
                        n = t ? t.length : 0;
                    for (this.clear(); ++e < n;) {
                        var r = t[e];
                        this.set(r[0], r[1])
                    }
                }

                function Ve() {
                    this.__data__ = []
                }

                function $e(t) {
                    var e = this.__data__,
                        n = gn(e, t);
                    return !(n < 0 || (n == e.length - 1 ? e.pop() : Vl.call(e, n, 1), 0))
                }

                function Ge(t) {
                    var e = this.__data__,
                        n = gn(e, t);
                    return n < 0 ? Q : e[n][1]
                }

                function Ke(t) {
                    return gn(this.__data__, t) > -1
                }

                function Xe(t, e) {
                    var n = this.__data__,
                        r = gn(n, t);
                    return r < 0 ? n.push([t, e]) : n[r][1] = e, this
                }

                function Ze(t) {
                    var e = -1,
                        n = t ? t.length : 0;
                    for (this.clear(); ++e < n;) {
                        var r = t[e];
                        this.set(r[0], r[1])
                    }
                }

                function Je() {
                    this.__data__ = {
                        hash: new Fe,
                        map: new(uc || Ue),
                        string: new Fe
                    }
                }

                function Ye(t) {
                    return qi(this, t)["delete"](t)
                }

                function Qe(t) {
                    return qi(this, t).get(t)
                }

                function tn(t) {
                    return qi(this, t).has(t)
                }

                function en(t, e) {
                    return qi(this, t).set(t, e), this
                }

                function nn(t) {
                    var e = -1,
                        n = t ? t.length : 0;
                    for (this.__data__ = new Ze; ++e < n;) this.add(t[e])
                }

                function rn(t) {
                    return this.__data__.set(t, rt), this
                }

                function on(t) {
                    return this.__data__.has(t)
                }

                function sn(t) {
                    this.__data__ = new Ue(t)
                }

                function un() {
                    this.__data__ = new Ue
                }

                function an(t) {
                    return this.__data__["delete"](t)
                }

                function ln(t) {
                    return this.__data__.get(t)
                }

                function cn(t) {
                    return this.__data__.has(t)
                }

                function fn(t, e) {
                    var n = this.__data__;
                    return n instanceof Ue && n.__data__.length == et && (n = this.__data__ = new Ze(n.__data__)), n.set(t, e), this
                }

                function hn(t, e, n, r) {
                    return t === Q || nu(t, jl[n]) && !Al.call(r, n) ? e : t
                }

                function dn(t, e, n) {
                    (n === Q || nu(t[e], n)) && ("number" != typeof e || n !== Q || e in t) || (t[e] = n)
                }

                function pn(t, e, n) {
                    var r = t[e];
                    Al.call(t, e) && nu(r, n) && (n !== Q || e in t) || (t[e] = n)
                }

                function gn(t, e) {
                    for (var n = t.length; n--;)
                        if (nu(t[n][0], e)) return n;
                    return -1
                }

                function vn(t, e, n, r) {
                    return Sc(t, function(t, i, o) {
                        e(r, t, n(t), o)
                    }), r
                }

                function mn(t, e) {
                    return t && ii(e, ia(e), t)
                }

                function yn(t, e) {
                    for (var n = -1, r = null == t, i = e.length, o = Array(i); ++n < i;) o[n] = r ? Q : ea(t, e[n]);
                    return o
                }

                function _n(t, e, n) {
                    return t === t && (n !== Q && (t = t <= n ? t : n), e !== Q && (t = t >= e ? t : e)), t
                }

                function wn(t, e, n, r, i, o, s) {
                    var u;
                    if (r && (u = o ? r(t, i, o, s) : r(t)), u !== Q) return u;
                    if (!yu(t)) return t;
                    var a = yf(t);
                    if (a) {
                        if (u = Vi(t), !e) return ri(t, u)
                    } else {
                        var l = Di(t),
                            f = l == Rt || l == qt;
                        if (_f(t)) return Vr(t, e);
                        if (l == Bt || l == Ct || f && !o) {
                            if (H(t)) return o ? t : {};
                            if (u = $i(f ? {} : t), !e) return oi(t, mn(u, t))
                        } else {
                            if (!zn[l]) return o ? t : {};
                            u = Gi(t, l, wn, e)
                        }
                    }
                    s || (s = new sn);
                    var h = s.get(t);
                    if (h) return h;
                    if (s.set(t, u), !a) var d = n ? Ti(t) : ia(t);
                    return c(d || t, function(i, o) {
                        d && (o = i, i = t[o]), pn(u, o, wn(i, e, n, r, o, t, s))
                    }), u
                }

                function An(t) {
                    var e = ia(t),
                        n = e.length;
                    return function(r) {
                        if (null == r) return !n;
                        for (var i = n; i--;) {
                            var o = e[i],
                                s = t[o],
                                u = r[o];
                            if (u === Q && !(o in Object(r)) || !s(u)) return !1
                        }
                        return !0
                    }
                }

                function Cn(t) {
                    return yu(t) ? Hl(t) : {}
                }

                function Tn(t, e, n) {
                    if ("function" != typeof t) throw new Sl(nt);
                    return $l(function() {
                        t.apply(Q, n)
                    }, e)
                }

                function In(t, e, n, r) {
                    var i = -1,
                        o = p,
                        s = !0,
                        u = t.length,
                        a = [],
                        l = e.length;
                    if (!u) return a;
                    n && (e = v(e, C(n))), r ? (o = g, s = !1) : e.length >= et && (o = I, s = !1, e = new nn(e));
                    t: for (; ++i < u;) {
                        var c = t[i],
                            f = n ? n(c) : c;
                        if (c = r || 0 !== c ? c : 0, s && f === f) {
                            for (var h = l; h--;)
                                if (e[h] === f) continue t;
                            a.push(c)
                        } else o(e, f, r) || a.push(c)
                    }
                    return a
                }

                function Rn(t, e) {
                    var n = !0;
                    return Sc(t, function(t, r, i) {
                        return n = !!e(t, r, i)
                    }), n
                }

                function qn(t, e, n) {
                    for (var r = -1, i = t.length; ++r < i;) {
                        var o = t[r],
                            s = e(o);
                        if (null != s && (u === Q ? s === s && !Iu(s) : n(s, u))) var u = s,
                            a = o
                    }
                    return a
                }

                function Pn(t, e, n, r) {
                    var i = t.length;
                    for (n = Bu(n), n < 0 && (n = -n > i ? 0 : i + n), r = r === Q || r > i ? i : Bu(r), r < 0 && (r += i), r = n > r ? 0 : Wu(r); n < r;) t[n++] = e;
                    return t
                }

                function Bn(t, e) {
                    var n = [];
                    return Sc(t, function(t, r, i) {
                        e(t, r, i) && n.push(t)
                    }), n
                }

                function Wn(t, e, n, r, i) {
                    var o = -1,
                        s = t.length;
                    for (n || (n = Xi), i || (i = []); ++o < s;) {
                        var u = t[o];
                        e > 0 && n(u) ? e > 1 ? Wn(u, e - 1, n, r, i) : m(i, u) : r || (i[i.length] = u)
                    }
                    return i
                }

                function Un(t, e) {
                    return t && jc(t, e, ia)
                }

                function Vn(t, e) {
                    return t && Oc(t, e, ia)
                }

                function $n(t, e) {
                    return d(e, function(e) {
                        return gu(t[e])
                    })
                }

                function Gn(t, e) {
                    e = Qi(e, t) ? [e] : Hr(e);
                    for (var n = 0, r = e.length; null != t && n < r;) t = t[co(e[n++])];
                    return n && n == r ? t : Q
                }

                function Kn(t, e, n) {
                    var r = e(t);
                    return yf(t) ? r : m(r, n(t))
                }

                function Xn(t, e) {
                    return t > e
                }

                function Zn(t, e) {
                    return null != t && (Al.call(t, e) || "object" == ("undefined" == typeof t ? "undefined" : o(t)) && e in t && null === Bi(t))
                }

                function Jn(t, e) {
                    return null != t && e in Object(t)
                }

                function Yn(t, e, n) {
                    return t >= tc(e, n) && t < Ql(e, n)
                }

                function Qn(t, e, n) {
                    for (var r = n ? g : p, i = t[0].length, o = t.length, s = o, u = Array(o), a = 1 / 0, l = []; s--;) {
                        var c = t[s];
                        s && e && (c = v(c, C(e))), a = tc(c.length, a), u[s] = !n && (e || i >= 120 && c.length >= 120) ? new nn(s && c) : Q
                    }
                    c = t[0];
                    var f = -1,
                        h = u[0];
                    t: for (; ++f < i && l.length < a;) {
                        var d = c[f],
                            m = e ? e(d) : d;
                        if (d = n || 0 !== d ? d : 0, !(h ? I(h, m) : r(l, m, n))) {
                            for (s = o; --s;) {
                                var y = u[s];
                                if (!(y ? I(y, m) : r(t[s], m, n))) continue t
                            }
                            h && h.push(m), l.push(d)
                        }
                    }
                    return l
                }

                function tr(t, e, n, r) {
                    return Un(t, function(t, i, o) {
                        e(r, n(t), i, o)
                    }), r
                }

                function er(t, e, n) {
                    Qi(e, t) || (e = Hr(e), t = ao(t, e), e = To(e));
                    var r = null == t ? t : t[co(e)];
                    return null == r ? Q : a(r, t, n)
                }

                function nr(t, e, n, r, i) {
                    return t === e || (null == t || null == e || !yu(t) && !_u(e) ? t !== t && e !== e : rr(t, e, nr, n, r, i))
                }

                function rr(t, e, n, r, i, o) {
                    var s = yf(t),
                        u = yf(e),
                        a = Tt,
                        l = Tt;
                    s || (a = Di(t), a = a == Ct ? Bt : a), u || (l = Di(e), l = l == Ct ? Bt : l);
                    var c = a == Bt && !H(t),
                        f = l == Bt && !H(e),
                        h = a == l;
                    if (h && !c) return o || (o = new sn), s || Nu(t) ? zi(t, e, n, r, i, o) : Ai(t, e, a, n, r, i, o);
                    if (!(i & vt)) {
                        var d = c && Al.call(t, "__wrapped__"),
                            p = f && Al.call(e, "__wrapped__");
                        if (d || p) {
                            var g = d ? t.value() : t,
                                v = p ? e.value() : e;
                            return o || (o = new sn), n(g, v, r, i, o)
                        }
                    }
                    return !!h && (o || (o = new sn), Ci(t, e, n, r, i, o))
                }

                function ir(t, e, n, r) {
                    var i = n.length,
                        o = i,
                        s = !r;
                    if (null == t) return !o;
                    for (t = Object(t); i--;) {
                        var u = n[i];
                        if (s && u[2] ? u[1] !== t[u[0]] : !(u[0] in t)) return !1
                    }
                    for (; ++i < o;) {
                        u = n[i];
                        var a = u[0],
                            l = t[a],
                            c = u[1];
                        if (s && u[2]) {
                            if (l === Q && !(a in t)) return !1
                        } else {
                            var f = new sn;
                            if (r) var h = r(l, c, a, t, e, f);
                            if (!(h === Q ? nr(c, l, r, gt | vt, f) : h)) return !1
                        }
                    }
                    return !0
                }

                function or(t) {
                    return !(!yu(t) || no(t)) && (gu(t) || H(t) ? Ll : Ae).test(fo(t))
                }

                function sr(t) {
                    return "function" == typeof t ? t : null == t ? Ga : "object" == ("undefined" == typeof t ? "undefined" : o(t)) ? yf(t) ? hr(t[0], t[1]) : fr(t) : el(t)
                }

                function ur(t) {
                    return Yl(Object(t))
                }

                function ar(t) {
                    t = null == t ? t : Object(t);
                    var e = [];
                    for (var n in t) e.push(n);
                    return e
                }

                function lr(t, e) {
                    return t < e
                }

                function cr(t, e) {
                    var n = -1,
                        r = ou(t) ? Array(t.length) : [];
                    return Sc(t, function(t, i, o) {
                        r[++n] = e(t, i, o)
                    }), r
                }

                function fr(t) {
                    var e = Fi(t);
                    return 1 == e.length && e[0][2] ? oo(e[0][0], e[0][1]) : function(n) {
                        return n === t || ir(n, t, e)
                    }
                }

                function hr(t, e) {
                    return Qi(t) && io(e) ? oo(co(t), e) : function(n) {
                        var r = ea(n, t);
                        return r === Q && r === e ? ra(n, t) : nr(e, r, Q, gt | vt)
                    }
                }

                function dr(t, e, n, r, i) {
                    if (t !== e) {
                        if (!yf(e) && !Nu(e)) var o = oa(e);
                        c(o || e, function(s, u) {
                            if (o && (u = s, s = e[u]), yu(s)) i || (i = new sn), pr(t, e, u, n, dr, r, i);
                            else {
                                var a = r ? r(t[u], s, u + "", t, e, i) : Q;
                                a === Q && (a = s), dn(t, u, a)
                            }
                        })
                    }
                }

                function pr(t, e, n, r, i, o, s) {
                    var u = t[n],
                        a = e[n],
                        l = s.get(a);
                    if (l) return void dn(t, n, l);
                    var c = o ? o(u, a, n + "", t, e, s) : Q,
                        f = c === Q;
                    f && (c = a, yf(a) || Nu(a) ? yf(u) ? c = u : su(u) ? c = ri(u) : (f = !1, c = wn(a, !0)) : Mu(a) || ru(a) ? ru(u) ? c = Hu(u) : !yu(u) || r && gu(u) ? (f = !1, c = wn(a, !0)) : c = u : f = !1), s.set(a, c), f && i(c, a, r, o, s), s["delete"](a), dn(t, n, c)
                }

                function gr(t, e) {
                    var n = t.length;
                    if (n) return e += e < 0 ? n : 0, Ji(e, n) ? t[e] : Q
                }

                function vr(t, e, n) {
                    var r = -1;
                    return e = v(e.length ? e : [Ga], C(Ri())), E(cr(t, function(t, n, i) {
                        return {
                            criteria: v(e, function(e) {
                                return e(t)
                            }),
                            index: ++r,
                            value: t
                        }
                    }), function(t, e) {
                        return ti(t, e, n)
                    })
                }

                function mr(t, e) {
                    return t = Object(t), y(e, function(e, n) {
                        return n in t && (e[n] = t[n]), e
                    }, {})
                }

                function yr(t, e) {
                    for (var n = -1, r = Ii(t), i = r.length, o = {}; ++n < i;) {
                        var s = r[n],
                            u = t[s];
                        e(u, s) && (o[s] = u)
                    }
                    return o
                }

                function _r(t) {
                    return function(e) {
                        return null == e ? Q : e[t]
                    }
                }

                function br(t) {
                    return function(e) {
                        return Gn(e, t)
                    }
                }

                function xr(t, e, n, r) {
                    var i = r ? k : S,
                        o = -1,
                        s = e.length,
                        u = t;
                    for (t === e && (e = ri(e)), n && (u = v(t, C(n))); ++o < s;)
                        for (var a = 0, l = e[o], c = n ? n(l) : l;
                            (a = i(u, c, a, r)) > -1;) u !== t && Vl.call(u, a, 1), Vl.call(t, a, 1);
                    return t
                }

                function wr(t, e) {
                    for (var n = t ? e.length : 0, r = n - 1; n--;) {
                        var i = e[n];
                        if (n == r || i !== o) {
                            var o = i;
                            if (Ji(i)) Vl.call(t, i, 1);
                            else if (Qi(i, t)) delete t[co(i)];
                            else {
                                var s = Hr(i),
                                    u = ao(t, s);
                                null != u && delete u[co(To(s))]
                            }
                        }
                    }
                    return t
                }

                function Sr(t, e) {
                    return t + Kl(nc() * (e - t + 1))
                }

                function kr(t, e, n, r) {
                    for (var i = -1, o = Ql(Gl((e - t) / (n || 1)), 0), s = Array(o); o--;) s[r ? o : ++i] = t, t += n;
                    return s
                }

                function jr(t, e) {
                    var n = "";
                    if (!t || e < 1 || e > jt) return n;
                    do e % 2 && (n += t), e = Kl(e / 2), e && (t += t); while (e);
                    return n
                }

                function Or(t, e, n, r) {
                    e = Qi(e, t) ? [e] : Hr(e);
                    for (var i = -1, o = e.length, s = o - 1, u = t; null != u && ++i < o;) {
                        var a = co(e[i]);
                        if (yu(u)) {
                            var l = n;
                            if (i != s) {
                                var c = u[a];
                                l = r ? r(c, a, u) : Q, l === Q && (l = null == c ? Ji(e[i + 1]) ? [] : {} : c)
                            }
                            pn(u, a, l)
                        }
                        u = u[a]
                    }
                    return t
                }

                function Er(t, e, n) {
                    var r = -1,
                        i = t.length;
                    e < 0 && (e = -e > i ? 0 : i + e), n = n > i ? i : n, n < 0 && (n += i), i = e > n ? 0 : n - e >>> 0, e >>>= 0;
                    for (var o = Array(i); ++r < i;) o[r] = t[r + e];
                    return o
                }

                function Mr(t, e) {
                    var n;
                    return Sc(t, function(t, r, i) {
                        return n = e(t, r, i), !n
                    }), !!n
                }

                function zr(t, e, n) {
                    var r = 0,
                        i = t ? t.length : r;
                    if ("number" == typeof e && e === e && i <= At) {
                        for (; r < i;) {
                            var o = r + i >>> 1,
                                s = t[o];
                            null !== s && !Iu(s) && (n ? s <= e : s < e) ? r = o + 1 : i = o
                        }
                        return i
                    }
                    return Ar(t, e, Ga, n)
                }

                function Ar(t, e, n, r) {
                    e = n(e);
                    for (var i = 0, o = t ? t.length : 0, s = e !== e, u = null === e, a = Iu(e), l = e === Q; i < o;) {
                        var c = Kl((i + o) / 2),
                            f = n(t[c]),
                            h = f !== Q,
                            d = null === f,
                            p = f === f,
                            g = Iu(f);
                        if (s) var v = r || p;
                        else v = l ? p && (r || h) : u ? p && h && (r || !d) : a ? p && h && !d && (r || !g) : !d && !g && (r ? f <= e : f < e);
                        v ? i = c + 1 : o = c
                    }
                    return tc(o, zt)
                }

                function Cr(t, e) {
                    for (var n = -1, r = t.length, i = 0, o = []; ++n < r;) {
                        var s = t[n],
                            u = e ? e(s) : s;
                        if (!n || !nu(u, a)) {
                            var a = u;
                            o[i++] = 0 === s ? 0 : s
                        }
                    }
                    return o
                }

                function Tr(t) {
                    return "number" == typeof t ? t : Iu(t) ? Et : +t
                }

                function Ir(t) {
                    if ("string" == typeof t) return t;
                    if (Iu(t)) return wc ? wc.call(t) : "";
                    var e = t + "";
                    return "0" == e && 1 / t == -kt ? "-0" : e
                }

                function Nr(t, e, n) {
                    var r = -1,
                        i = p,
                        o = t.length,
                        s = !0,
                        u = [],
                        a = u;
                    if (n) s = !1, i = g;
                    else if (o >= et) {
                        var l = e ? null : Mc(t);
                        if (l) return G(l);
                        s = !1, i = I, a = new nn
                    } else a = e ? [] : u;
                    t: for (; ++r < o;) {
                        var c = t[r],
                            f = e ? e(c) : c;
                        if (c = n || 0 !== c ? c : 0, s && f === f) {
                            for (var h = a.length; h--;)
                                if (a[h] === f) continue t;
                            e && a.push(f), u.push(c)
                        } else i(a, f, n) || (a !== u && a.push(f), u.push(c))
                    }
                    return u
                }

                function Lr(t, e) {
                    e = Qi(e, t) ? [e] : Hr(e), t = ao(t, e);
                    var n = co(To(e));
                    return !(null != t && Zn(t, n)) || delete t[n]
                }

                function Rr(t, e, n, r) {
                    return Or(t, e, n(Gn(t, e)), r)
                }

                function qr(t, e, n, r) {
                    for (var i = t.length, o = r ? i : -1;
                        (r ? o-- : ++o < i) && e(t[o], o, t););
                    return n ? Er(t, r ? 0 : o, r ? o + 1 : i) : Er(t, r ? o + 1 : 0, r ? i : o)
                }

                function Fr(t, e) {
                    var n = t;
                    return n instanceof i && (n = n.value()), y(e, function(t, e) {
                        return e.func.apply(e.thisArg, m([t], e.args))
                    }, n)
                }

                function Pr(t, e, n) {
                    for (var r = -1, i = t.length; ++r < i;) var o = o ? m(In(o, t[r], e, n), In(t[r], o, e, n)) : t[r];
                    return o && o.length ? Nr(o, e, n) : []
                }

                function Br(t, e, n) {
                    for (var r = -1, i = t.length, o = e.length, s = {}; ++r < i;) {
                        var u = r < o ? e[r] : Q;
                        n(s, t[r], u)
                    }
                    return s
                }

                function Wr(t) {
                    return su(t) ? t : []
                }

                function Dr(t) {
                    return "function" == typeof t ? t : Ga
                }

                function Hr(t) {
                    return yf(t) ? t : Nc(t)
                }

                function Ur(t, e, n) {
                    var r = t.length;
                    return n = n === Q ? r : n, !e && n >= r ? t : Er(t, e, n)
                }

                function Vr(t, e) {
                    if (e) return t.slice();
                    var n = new t.constructor(t.length);
                    return t.copy(n), n
                }

                function $r(t) {
                    var e = new t.constructor(t.byteLength);
                    return new Pl(e).set(new Pl(t)), e
                }

                function Gr(t, e) {
                    var n = e ? $r(t.buffer) : t.buffer;
                    return new t.constructor(n, t.byteOffset, t.byteLength)
                }

                function Kr(t, e, n) {
                    return y(e ? n(V(t), !0) : V(t), s, new t.constructor)
                }

                function Xr(t) {
                    var e = new t.constructor(t.source, Oe.exec(t));
                    return e.lastIndex = t.lastIndex, e
                }

                function Zr(t, e, n) {
                    return y(e ? n(G(t), !0) : G(t), u, new t.constructor)
                }

                function Jr(t) {
                    return xc ? Object(xc.call(t)) : {}
                }

                function Yr(t, e) {
                    var n = e ? $r(t.buffer) : t.buffer;
                    return new t.constructor(n, t.byteOffset, t.length)
                }

                function Qr(t, e) {
                    if (t !== e) {
                        var n = t !== Q,
                            r = null === t,
                            i = t === t,
                            o = Iu(t),
                            s = e !== Q,
                            u = null === e,
                            a = e === e,
                            l = Iu(e);
                        if (!u && !l && !o && t > e || o && s && a && !u && !l || r && s && a || !n && a || !i) return 1;
                        if (!r && !o && !l && t < e || l && n && i && !r && !o || u && n && i || !s && i || !a) return -1
                    }
                    return 0
                }

                function ti(t, e, n) {
                    for (var r = -1, i = t.criteria, o = e.criteria, s = i.length, u = n.length; ++r < s;) {
                        var a = Qr(i[r], o[r]);
                        if (a) return r >= u ? a : a * ("desc" == n[r] ? -1 : 1)
                    }
                    return t.index - e.index
                }

                function ei(t, e, n, r) {
                    for (var i = -1, o = t.length, s = n.length, u = -1, a = e.length, l = Ql(o - s, 0), c = Array(a + l), f = !r; ++u < a;) c[u] = e[u];
                    for (; ++i < s;)(f || i < o) && (c[n[i]] = t[i]);
                    for (; l--;) c[u++] = t[i++];
                    return c
                }

                function ni(t, e, n, r) {
                    for (var i = -1, o = t.length, s = -1, u = n.length, a = -1, l = e.length, c = Ql(o - u, 0), f = Array(c + l), h = !r; ++i < c;) f[i] = t[i];
                    for (var d = i; ++a < l;) f[d + a] = e[a];
                    for (; ++s < u;)(h || i < o) && (f[d + n[s]] = t[i++]);
                    return f
                }

                function ri(t, e) {
                    var n = -1,
                        r = t.length;
                    for (e || (e = Array(r)); ++n < r;) e[n] = t[n];
                    return e
                }

                function ii(t, e, n, r) {
                    n || (n = {});
                    for (var i = -1, o = e.length; ++i < o;) {
                        var s = e[i];
                        pn(n, s, r ? r(n[s], t[s], s, n, t) : t[s])
                    }
                    return n
                }

                function oi(t, e) {
                    return ii(t, Wi(t), e)
                }

                function si(t, e) {
                    return function(n, r) {
                        var i = yf(n) ? l : vn,
                            o = e ? e() : {};
                        return i(n, t, Ri(r), o)
                    }
                }

                function ui(t) {
                    return $s(function(e, n) {
                        var r = -1,
                            i = n.length,
                            o = i > 1 ? n[i - 1] : Q,
                            s = i > 2 ? n[2] : Q;
                        for (o = t.length > 3 && "function" == typeof o ? (i--, o) : Q, s && Yi(n[0], n[1], s) && (o = i < 3 ? Q : o, i = 1), e = Object(e); ++r < i;) {
                            var u = n[r];
                            u && t(e, u, r, o)
                        }
                        return e
                    })
                }

                function ai(t, e) {
                    return function(n, r) {
                        if (null == n) return n;
                        if (!ou(n)) return t(n, r);
                        for (var i = n.length, o = e ? i : -1, s = Object(n);
                            (e ? o-- : ++o < i) && r(s[o], o, s) !== !1;);
                        return n
                    }
                }

                function li(t) {
                    return function(e, n, r) {
                        for (var i = -1, o = Object(e), s = r(e), u = s.length; u--;) {
                            var a = s[t ? u : ++i];
                            if (n(o[a], a, o) === !1) break
                        }
                        return e
                    }
                }

                function ci(t, e, n) {
                    function r() {
                        return (this && this !== Dn && this instanceof r ? o : t).apply(i ? n : this, arguments)
                    }
                    var i = e & ot,
                        o = di(t);
                    return r
                }

                function fi(t) {
                    return function(e) {
                        e = Vu(e);
                        var n = kn.test(e) ? Z(e) : Q,
                            r = n ? n[0] : e.charAt(0),
                            i = n ? Ur(n, 1).join("") : e.slice(1);
                        return r[t]() + i
                    }
                }

                function hi(t) {
                    return function(e) {
                        return y(Ha(Sa(e).replace(bn, "")), t, "")
                    }
                }

                function di(t) {
                    return function() {
                        var e = arguments;
                        switch (e.length) {
                            case 0:
                                return new t;
                            case 1:
                                return new t(e[0]);
                            case 2:
                                return new t(e[0], e[1]);
                            case 3:
                                return new t(e[0], e[1], e[2]);
                            case 4:
                                return new t(e[0], e[1], e[2], e[3]);
                            case 5:
                                return new t(e[0], e[1], e[2], e[3], e[4]);
                            case 6:
                                return new t(e[0], e[1], e[2], e[3], e[4], e[5]);
                            case 7:
                                return new t(e[0], e[1], e[2], e[3], e[4], e[5], e[6])
                        }
                        var n = Cn(t.prototype),
                            r = t.apply(n, e);
                        return yu(r) ? r : n
                    }
                }

                function pi(t, e, n) {
                    function r() {
                        for (var o = arguments.length, s = Array(o), u = o, l = Li(r); u--;) s[u] = arguments[u];
                        var c = o < 3 && s[0] !== l && s[o - 1] !== l ? [] : $(s, l);
                        return o -= c.length, o < n ? ji(t, e, mi, r.placeholder, Q, s, c, Q, Q, n - o) : a(this && this !== Dn && this instanceof r ? i : t, this, s)
                    }
                    var i = di(t);
                    return r
                }

                function gi(t) {
                    return function(e, n, r) {
                        var i = Object(e);
                        if (n = Ri(n, 3), !ou(e)) var o = ia(e);
                        var s = t(o || e, function(t, e) {
                            return o && (e = t, t = i[e]), n(t, e, i)
                        }, r);
                        return s > -1 ? e[o ? o[s] : s] : Q
                    }
                }

                function vi(t) {
                    return $s(function(e) {
                        e = Wn(e, 1);
                        var n = e.length,
                            i = n,
                            o = r.prototype.thru;
                        for (t && e.reverse(); i--;) {
                            var s = e[i];
                            if ("function" != typeof s) throw new Sl(nt);
                            if (o && !u && "wrapper" == Ni(s)) var u = new r([], (!0))
                        }
                        for (i = u ? i : n; ++i < n;) {
                            s = e[i];
                            var a = Ni(s),
                                l = "wrapper" == a ? zc(s) : Q;
                            u = l && eo(l[0]) && l[1] == (ht | at | ct | dt) && !l[4].length && 1 == l[9] ? u[Ni(l[0])].apply(u, l[3]) : 1 == s.length && eo(s) ? u[a]() : u.thru(s)
                        }
                        return function() {
                            var t = arguments,
                                r = t[0];
                            if (u && 1 == t.length && yf(r) && r.length >= et) return u.plant(r).value();
                            for (var i = 0, o = n ? e[i].apply(this, t) : r; ++i < n;) o = e[i].call(this, o);
                            return o
                        }
                    })
                }

                function mi(t, e, n, r, i, o, s, u, a, l) {
                    function c() {
                        for (var m = arguments.length, y = Array(m), _ = m; _--;) y[_] = arguments[_];
                        if (p) var b = Li(c),
                            x = q(y, b);
                        if (r && (y = ei(y, r, i, p)), o && (y = ni(y, o, s, p)), m -= x, p && m < l) {
                            var w = $(y, b);
                            return ji(t, e, mi, c.placeholder, n, y, w, u, a, l - m)
                        }
                        var S = h ? n : this,
                            k = d ? S[t] : t;
                        return m = y.length, u ? y = lo(y, u) : g && m > 1 && y.reverse(), f && a < m && (y.length = a), this && this !== Dn && this instanceof c && (k = v || di(k)), k.apply(S, y)
                    }
                    var f = e & ht,
                        h = e & ot,
                        d = e & st,
                        p = e & (at | lt),
                        g = e & pt,
                        v = d ? Q : di(t);
                    return c
                }

                function yi(t, e) {
                    return function(n, r) {
                        return tr(n, t, e(r), {})
                    }
                }

                function _i(t) {
                    return function(e, n) {
                        var r;
                        if (e === Q && n === Q) return 0;
                        if (e !== Q && (r = e), n !== Q) {
                            if (r === Q) return n;
                            "string" == typeof e || "string" == typeof n ? (e = Ir(e), n = Ir(n)) : (e = Tr(e), n = Tr(n)), r = t(e, n)
                        }
                        return r
                    }
                }

                function bi(t) {
                    return $s(function(e) {
                        return e = 1 == e.length && yf(e[0]) ? v(e[0], C(Ri())) : v(Wn(e, 1, Zi), C(Ri())), $s(function(n) {
                            var r = this;
                            return t(e, function(t) {
                                return a(t, r, n)
                            })
                        })
                    })
                }

                function xi(t, e) {
                    e = e === Q ? " " : Ir(e);
                    var n = e.length;
                    if (n < 2) return n ? jr(e, t) : e;
                    var r = jr(e, Gl(t / X(e)));
                    return kn.test(e) ? Ur(Z(r), 0, t).join("") : r.slice(0, t)
                }

                function wi(t, e, n, r) {
                    function i() {
                        for (var e = -1, u = arguments.length, l = -1, c = r.length, f = Array(c + u), h = this && this !== Dn && this instanceof i ? s : t; ++l < c;) f[l] = r[l];
                        for (; u--;) f[l++] = arguments[++e];
                        return a(h, o ? n : this, f)
                    }
                    var o = e & ot,
                        s = di(t);
                    return i
                }

                function Si(t) {
                    return function(e, n, r) {
                        return r && "number" != typeof r && Yi(e, n, r) && (n = r = Q), e = Du(e), e = e === e ? e : 0, n === Q ? (n = e, e = 0) : n = Du(n) || 0, r = r === Q ? e < n ? 1 : -1 : Du(r) || 0, kr(e, n, r, t)
                    }
                }

                function ki(t) {
                    return function(e, n) {
                        return "string" == typeof e && "string" == typeof n || (e = Du(e), n = Du(n)), t(e, n)
                    }
                }

                function ji(t, e, n, r, i, o, s, u, a, l) {
                    var c = e & at,
                        f = c ? s : Q,
                        h = c ? Q : s,
                        d = c ? o : Q,
                        p = c ? Q : o;
                    e |= c ? ct : ft, e &= ~(c ? ft : ct), e & ut || (e &= ~(ot | st));
                    var g = [t, e, i, d, f, p, h, u, a, l],
                        v = n.apply(Q, g);
                    return eo(t) && Ic(v, g), v.placeholder = r, v
                }

                function Oi(t) {
                    var e = xl[t];
                    return function(t, n) {
                        if (t = Du(t), n = tc(Bu(n), 292)) {
                            var r = (Vu(t) + "e").split("e");
                            return r = (Vu(e(r[0] + "e" + (+r[1] + n))) + "e").split("e"), +(r[0] + "e" + (+r[1] - n))
                        }
                        return e(t)
                    }
                }

                function Ei(t) {
                    return function(e) {
                        var n = Di(e);
                        return n == Ft ? V(e) : n == Ht ? K(e) : A(e, t(e))
                    }
                }

                function Mi(t, e, n, r, i, o, s, u) {
                    var a = e & st;
                    if (!a && "function" != typeof t) throw new Sl(nt);
                    var l = r ? r.length : 0;
                    if (l || (e &= ~(ct | ft), r = i = Q), s = s === Q ? s : Ql(Bu(s), 0), u = u === Q ? u : Bu(u), l -= i ? i.length : 0, e & ft) {
                        var c = r,
                            f = i;
                        r = i = Q
                    }
                    var h = a ? Q : zc(t),
                        d = [t, e, n, r, i, c, f, o, s, u];
                    if (h && so(d, h), t = d[0], e = d[1], n = d[2], r = d[3], i = d[4], u = d[9] = null == d[9] ? a ? 0 : t.length : Ql(d[9] - l, 0), !u && e & (at | lt) && (e &= ~(at | lt)), e && e != ot) p = e == at || e == lt ? pi(t, e, u) : e != ct && e != (ot | ct) || i.length ? mi.apply(Q, d) : wi(t, e, n, r);
                    else var p = ci(t, e, n);
                    return (h ? Ec : Ic)(p, d)
                }

                function zi(t, e, n, r, i, o) {
                    var s = i & vt,
                        u = t.length,
                        a = e.length;
                    if (u != a && !(s && a > u)) return !1;
                    var l = o.get(t);
                    if (l) return l == e;
                    var c = -1,
                        f = !0,
                        h = i & gt ? new nn : Q;
                    for (o.set(t, e); ++c < u;) {
                        var d = t[c],
                            p = e[c];
                        if (r) var g = s ? r(p, d, c, e, t, o) : r(d, p, c, t, e, o);
                        if (g !== Q) {
                            if (g) continue;
                            f = !1;
                            break
                        }
                        if (h) {
                            if (!b(e, function(t, e) {
                                    if (!h.has(e) && (d === t || n(d, t, r, i, o))) return h.add(e)
                                })) {
                                f = !1;
                                break
                            }
                        } else if (d !== p && !n(d, p, r, i, o)) {
                            f = !1;
                            break
                        }
                    }
                    return o["delete"](t), f
                }

                function Ai(t, e, n, r, i, o, s) {
                    switch (n) {
                        case Xt:
                            if (t.byteLength != e.byteLength || t.byteOffset != e.byteOffset) return !1;
                            t = t.buffer, e = e.buffer;
                        case Kt:
                            return !(t.byteLength != e.byteLength || !r(new Pl(t), new Pl(e)));
                        case It:
                        case Nt:
                            return +t == +e;
                        case Lt:
                            return t.name == e.name && t.message == e.message;
                        case Pt:
                            return t != +t ? e != +e : t == +e;
                        case Dt:
                        case Ut:
                            return t == e + "";
                        case Ft:
                            var u = V;
                        case Ht:
                            var a = o & vt;
                            if (u || (u = G), t.size != e.size && !a) return !1;
                            var l = s.get(t);
                            return l ? l == e : (o |= gt, s.set(t, e), zi(u(t), u(e), r, i, o, s));
                        case Vt:
                            if (xc) return xc.call(t) == xc.call(e)
                    }
                    return !1
                }

                function Ci(t, e, n, r, i, o) {
                    var s = i & vt,
                        u = ia(t),
                        a = u.length;
                    if (a != ia(e).length && !s) return !1;
                    for (var l = a; l--;) {
                        var c = u[l];
                        if (!(s ? c in e : Zn(e, c))) return !1
                    }
                    var f = o.get(t);
                    if (f) return f == e;
                    var h = !0;
                    o.set(t, e);
                    for (var d = s; ++l < a;) {
                        c = u[l];
                        var p = t[c],
                            g = e[c];
                        if (r) var v = s ? r(g, p, c, e, t, o) : r(p, g, c, t, e, o);
                        if (!(v === Q ? p === g || n(p, g, r, i, o) : v)) {
                            h = !1;
                            break
                        }
                        d || (d = "constructor" == c)
                    }
                    if (h && !d) {
                        var m = t.constructor,
                            y = e.constructor;
                        m != y && "constructor" in t && "constructor" in e && !("function" == typeof m && m instanceof m && "function" == typeof y && y instanceof y) && (h = !1)
                    }
                    return o["delete"](t), h
                }

                function Ti(t) {
                    return Kn(t, ia, Wi)
                }

                function Ii(t) {
                    return Kn(t, oa, Cc)
                }

                function Ni(t) {
                    for (var e = t.name + "", n = pc[e], r = Al.call(pc, e) ? n.length : 0; r--;) {
                        var i = n[r],
                            o = i.func;
                        if (null == o || o == t) return i.name
                    }
                    return e
                }

                function Li(t) {
                    return (Al.call(e, "placeholder") ? e : t).placeholder
                }

                function Ri() {
                    var t = e.iteratee || Ka;
                    return t = t === Ka ? sr : t, arguments.length ? t(arguments[0], arguments[1]) : t
                }

                function qi(t, e) {
                    var n = t.__data__;
                    return to(e) ? n["string" == typeof e ? "string" : "hash"] : n.map
                }

                function Fi(t) {
                    for (var e = ia(t), n = e.length; n--;) {
                        var r = e[n],
                            i = t[r];
                        e[n] = [r, i, io(i)]
                    }
                    return e
                }

                function Pi(t, e) {
                    var n = W(t, e);
                    return or(n) ? n : Q
                }

                function Bi(t) {
                    return Xl(Object(t))
                }

                function Wi(t) {
                    return Wl(Object(t))
                }

                function Di(t) {
                    return Il.call(t)
                }

                function Hi(t, e, n) {
                    for (var r = -1, i = n.length; ++r < i;) {
                        var o = n[r],
                            s = o.size;
                        switch (o.type) {
                            case "drop":
                                t += s;
                                break;
                            case "dropRight":
                                e -= s;
                                break;
                            case "take":
                                e = tc(e, t + s);
                                break;
                            case "takeRight":
                                t = Ql(t, e - s)
                        }
                    }
                    return {
                        start: t,
                        end: e
                    }
                }

                function Ui(t, e, n) {
                    e = Qi(e, t) ? [e] : Hr(e);
                    for (var r, i = -1, o = e.length; ++i < o;) {
                        var s = co(e[i]);
                        if (!(r = null != t && n(t, s))) break;
                        t = t[s]
                    }
                    if (r) return r;
                    var o = t ? t.length : 0;
                    return !!o && mu(o) && Ji(s, o) && (yf(t) || Tu(t) || ru(t))
                }

                function Vi(t) {
                    var e = t.length,
                        n = t.constructor(e);
                    return e && "string" == typeof t[0] && Al.call(t, "index") && (n.index = t.index, n.input = t.input), n
                }

                function $i(t) {
                    return "function" != typeof t.constructor || ro(t) ? {} : Cn(Bi(t))
                }

                function Gi(t, e, n, r) {
                    var i = t.constructor;
                    switch (e) {
                        case Kt:
                            return $r(t);
                        case It:
                        case Nt:
                            return new i((+t));
                        case Xt:
                            return Gr(t, r);
                        case Zt:
                        case Jt:
                        case Yt:
                        case Qt:
                        case te:
                        case ee:
                        case ne:
                        case re:
                        case ie:
                            return Yr(t, r);
                        case Ft:
                            return Kr(t, r, n);
                        case Pt:
                        case Ut:
                            return new i(t);
                        case Dt:
                            return Xr(t);
                        case Ht:
                            return Zr(t, r, n);
                        case Vt:
                            return Jr(t)
                    }
                }

                function Ki(t) {
                    var e = t ? t.length : Q;
                    return mu(e) && (yf(t) || Tu(t) || ru(t)) ? z(e, String) : null
                }

                function Xi(t) {
                    return yf(t) || ru(t)
                }

                function Zi(t) {
                    return yf(t) && !(2 == t.length && !gu(t[0]))
                }

                function Ji(t, e) {
                    return e = null == e ? jt : e, !!e && ("number" == typeof t || Te.test(t)) && t > -1 && t % 1 == 0 && t < e
                }

                function Yi(t, e, n) {
                    if (!yu(n)) return !1;
                    var r = "undefined" == typeof e ? "undefined" : o(e);
                    return !!("number" == r ? ou(n) && Ji(e, n.length) : "string" == r && e in n) && nu(n[e], t)
                }

                function Qi(t, e) {
                    if (yf(t)) return !1;
                    var n = "undefined" == typeof t ? "undefined" : o(t);
                    return !("number" != n && "symbol" != n && "boolean" != n && null != t && !Iu(t)) || ve.test(t) || !ge.test(t) || null != e && t in Object(e)
                }

                function to(t) {
                    var e = "undefined" == typeof t ? "undefined" : o(t);
                    return "string" == e || "number" == e || "symbol" == e || "boolean" == e ? "__proto__" !== t : null === t
                }

                function eo(t) {
                    var n = Ni(t),
                        r = e[n];
                    if ("function" != typeof r || !(n in i.prototype)) return !1;
                    if (t === r) return !0;
                    var o = zc(r);
                    return !!o && t === o[0]
                }

                function no(t) {
                    return !!Ml && Ml in t
                }

                function ro(t) {
                    var e = t && t.constructor;
                    return t === ("function" == typeof e && e.prototype || jl)
                }

                function io(t) {
                    return t === t && !yu(t)
                }

                function oo(t, e) {
                    return function(n) {
                        return null != n && n[t] === e && (e !== Q || t in Object(n))
                    }
                }

                function so(t, e) {
                    var n = t[1],
                        r = e[1],
                        i = n | r,
                        o = i < (ot | st | ht),
                        s = r == ht && n == at || r == ht && n == dt && t[7].length <= e[8] || r == (ht | dt) && e[7].length <= e[8] && n == at;
                    if (!o && !s) return t;
                    r & ot && (t[2] = e[2], i |= n & ot ? 0 : ut);
                    var u = e[3];
                    if (u) {
                        var a = t[3];
                        t[3] = a ? ei(a, u, e[4]) : u, t[4] = a ? $(t[3], it) : e[4]
                    }
                    return u = e[5], u && (a = t[5], t[5] = a ? ni(a, u, e[6]) : u, t[6] = a ? $(t[5], it) : e[6]), u = e[7], u && (t[7] = u), r & ht && (t[8] = null == t[8] ? e[8] : tc(t[8], e[8])), null == t[9] && (t[9] = e[9]), t[0] = e[0], t[1] = i, t
                }

                function uo(t, e, n, r, i, o) {
                    return yu(t) && yu(e) && dr(t, e, Q, uo, o.set(e, t)), t
                }

                function ao(t, e) {
                    return 1 == e.length ? t : Gn(t, Er(e, 0, -1))
                }

                function lo(t, e) {
                    for (var n = t.length, r = tc(e.length, n), i = ri(t); r--;) {
                        var o = e[r];
                        t[r] = Ji(o, n) ? i[o] : Q
                    }
                    return t
                }

                function co(t) {
                    if ("string" == typeof t || Iu(t)) return t;
                    var e = t + "";
                    return "0" == e && 1 / t == -kt ? "-0" : e
                }

                function fo(t) {
                    if (null != t) {
                        try {
                            return zl.call(t)
                        } catch (e) {}
                        try {
                            return t + ""
                        } catch (e) {}
                    }
                    return ""
                }

                function ho(t) {
                    if (t instanceof i) return t.clone();
                    var e = new r(t.__wrapped__, t.__chain__);
                    return e.__actions__ = ri(t.__actions__), e.__index__ = t.__index__, e.__values__ = t.__values__, e
                }

                function po(t, e, n) {
                    e = (n ? Yi(t, e, n) : e === Q) ? 1 : Ql(Bu(e), 0);
                    var r = t ? t.length : 0;
                    if (!r || e < 1) return [];
                    for (var i = 0, o = 0, s = Array(Gl(r / e)); i < r;) s[o++] = Er(t, i, i += e);
                    return s
                }

                function go(t) {
                    for (var e = -1, n = t ? t.length : 0, r = 0, i = []; ++e < n;) {
                        var o = t[e];
                        o && (i[r++] = o)
                    }
                    return i
                }

                function vo() {
                    for (var t = arguments.length, e = Array(t ? t - 1 : 0), n = arguments[0], r = t; r--;) e[r - 1] = arguments[r];
                    return t ? m(yf(n) ? ri(n) : [n], Wn(e, 1)) : []
                }

                function mo(t, e, n) {
                    var r = t ? t.length : 0;
                    return r ? (e = n || e === Q ? 1 : Bu(e), Er(t, e < 0 ? 0 : e, r)) : []
                }

                function yo(t, e, n) {
                    var r = t ? t.length : 0;
                    return r ? (e = n || e === Q ? 1 : Bu(e), e = r - e, Er(t, 0, e < 0 ? 0 : e)) : []
                }

                function _o(t, e) {
                    return t && t.length ? qr(t, Ri(e, 3), !0, !0) : []
                }

                function bo(t, e) {
                    return t && t.length ? qr(t, Ri(e, 3), !0) : []
                }

                function xo(t, e, n, r) {
                    var i = t ? t.length : 0;
                    return i ? (n && "number" != typeof n && Yi(t, e, n) && (n = 0, r = i), Pn(t, e, n, r)) : []
                }

                function wo(t, e, n) {
                    var r = t ? t.length : 0;
                    if (!r) return -1;
                    var i = null == n ? 0 : Bu(n);
                    return i < 0 && (i = Ql(r + i, 0)), w(t, Ri(e, 3), i)
                }

                function So(t, e, n) {
                    var r = t ? t.length : 0;
                    if (!r) return -1;
                    var i = r - 1;
                    return n !== Q && (i = Bu(n), i = n < 0 ? Ql(r + i, 0) : tc(i, r - 1)), w(t, Ri(e, 3), i, !0)
                }

                function ko(t) {
                    return (t ? t.length : 0) ? Wn(t, 1) : []
                }

                function jo(t) {
                    return (t ? t.length : 0) ? Wn(t, kt) : []
                }

                function Oo(t, e) {
                    return (t ? t.length : 0) ? (e = e === Q ? 1 : Bu(e), Wn(t, e)) : []
                }

                function Eo(t) {
                    for (var e = -1, n = t ? t.length : 0, r = {}; ++e < n;) {
                        var i = t[e];
                        r[i[0]] = i[1]
                    }
                    return r
                }

                function Mo(t) {
                    return t && t.length ? t[0] : Q
                }

                function zo(t, e, n) {
                    var r = t ? t.length : 0;
                    if (!r) return -1;
                    var i = null == n ? 0 : Bu(n);
                    return i < 0 && (i = Ql(r + i, 0)), S(t, e, i)
                }

                function Ao(t) {
                    return yo(t, 1)
                }

                function Co(t, e) {
                    return t ? Jl.call(t, e) : ""
                }

                function To(t) {
                    var e = t ? t.length : 0;
                    return e ? t[e - 1] : Q
                }

                function Io(t, e, n) {
                    var r = t ? t.length : 0;
                    if (!r) return -1;
                    var i = r;
                    if (n !== Q && (i = Bu(n), i = (i < 0 ? Ql(r + i, 0) : tc(i, r - 1)) + 1), e !== e) return D(t, i - 1, !0);
                    for (; i--;)
                        if (t[i] === e) return i;
                    return -1
                }

                function No(t, e) {
                    return t && t.length ? gr(t, Bu(e)) : Q
                }

                function Lo(t, e) {
                    return t && t.length && e && e.length ? xr(t, e) : t
                }

                function Ro(t, e, n) {
                    return t && t.length && e && e.length ? xr(t, e, Ri(n)) : t
                }

                function qo(t, e, n) {
                    return t && t.length && e && e.length ? xr(t, e, Q, n) : t
                }

                function Fo(t, e) {
                    var n = [];
                    if (!t || !t.length) return n;
                    var r = -1,
                        i = [],
                        o = t.length;
                    for (e = Ri(e, 3); ++r < o;) {
                        var s = t[r];
                        e(s, r, t) && (n.push(s), i.push(r))
                    }
                    return wr(t, i), n
                }

                function Po(t) {
                    return t ? ic.call(t) : t
                }

                function Bo(t, e, n) {
                    var r = t ? t.length : 0;
                    return r ? (n && "number" != typeof n && Yi(t, e, n) ? (e = 0, n = r) : (e = null == e ? 0 : Bu(e), n = n === Q ? r : Bu(n)), Er(t, e, n)) : []
                }

                function Wo(t, e) {
                    return zr(t, e)
                }

                function Do(t, e, n) {
                    return Ar(t, e, Ri(n))
                }

                function Ho(t, e) {
                    var n = t ? t.length : 0;
                    if (n) {
                        var r = zr(t, e);
                        if (r < n && nu(t[r], e)) return r
                    }
                    return -1
                }

                function Uo(t, e) {
                    return zr(t, e, !0)
                }

                function Vo(t, e, n) {
                    return Ar(t, e, Ri(n), !0)
                }

                function $o(t, e) {
                    if (t ? t.length : 0) {
                        var n = zr(t, e, !0) - 1;
                        if (nu(t[n], e)) return n
                    }
                    return -1
                }

                function Go(t) {
                    return t && t.length ? Cr(t) : []
                }

                function Ko(t, e) {
                    return t && t.length ? Cr(t, Ri(e)) : []
                }

                function Xo(t) {
                    return mo(t, 1)
                }

                function Zo(t, e, n) {
                    return t && t.length ? (e = n || e === Q ? 1 : Bu(e), Er(t, 0, e < 0 ? 0 : e)) : []
                }

                function Jo(t, e, n) {
                    var r = t ? t.length : 0;
                    return r ? (e = n || e === Q ? 1 : Bu(e), e = r - e, Er(t, e < 0 ? 0 : e, r)) : []
                }

                function Yo(t, e) {
                    return t && t.length ? qr(t, Ri(e, 3), !1, !0) : []
                }

                function Qo(t, e) {
                    return t && t.length ? qr(t, Ri(e, 3)) : []
                }

                function ts(t) {
                    return t && t.length ? Nr(t) : []
                }

                function es(t, e) {
                    return t && t.length ? Nr(t, Ri(e)) : []
                }

                function ns(t, e) {
                    return t && t.length ? Nr(t, Q, e) : []
                }

                function rs(t) {
                    if (!t || !t.length) return [];
                    var e = 0;
                    return t = d(t, function(t) {
                        if (su(t)) return e = Ql(t.length, e), !0
                    }), z(e, function(e) {
                        return v(t, _r(e))
                    })
                }

                function is(t, e) {
                    if (!t || !t.length) return [];
                    var n = rs(t);
                    return null == e ? n : v(n, function(t) {
                        return a(e, Q, t)
                    })
                }

                function os(t, e) {
                    return Br(t || [], e || [], pn)
                }

                function ss(t, e) {
                    return Br(t || [], e || [], Or)
                }

                function us(t) {
                    var n = e(t);
                    return n.__chain__ = !0, n
                }

                function as(t, e) {
                    return e(t), t
                }

                function ls(t, e) {
                    return e(t)
                }

                function cs() {
                    return us(this)
                }

                function fs() {
                    return new r(this.value(), this.__chain__)
                }

                function hs() {
                    this.__values__ === Q && (this.__values__ = Fu(this.value()));
                    var t = this.__index__ >= this.__values__.length;
                    return {
                        done: t,
                        value: t ? Q : this.__values__[this.__index__++]
                    }
                }

                function ds() {
                    return this
                }

                function ps(t) {
                    for (var e, r = this; r instanceof n;) {
                        var i = ho(r);
                        i.__index__ = 0, i.__values__ = Q, e ? o.__wrapped__ = i : e = i;
                        var o = i;
                        r = r.__wrapped__
                    }
                    return o.__wrapped__ = t, e
                }

                function gs() {
                    var t = this.__wrapped__;
                    if (t instanceof i) {
                        var e = t;
                        return this.__actions__.length && (e = new i(this)), e = e.reverse(), e.__actions__.push({
                            func: ls,
                            args: [Po],
                            thisArg: Q
                        }), new r(e, this.__chain__)
                    }
                    return this.thru(Po)
                }

                function vs() {
                    return Fr(this.__wrapped__, this.__actions__)
                }

                function ms(t, e, n) {
                    var r = yf(t) ? h : Rn;
                    return n && Yi(t, e, n) && (e = Q), r(t, Ri(e, 3))
                }

                function ys(t, e) {
                    return (yf(t) ? d : Bn)(t, Ri(e, 3))
                }

                function _s(t, e) {
                    return Wn(js(t, e), 1)
                }

                function bs(t, e) {
                    return Wn(js(t, e), kt)
                }

                function xs(t, e, n) {
                    return n = n === Q ? 1 : Bu(n), Wn(js(t, e), n)
                }

                function ws(t, e) {
                    return (yf(t) ? c : Sc)(t, Ri(e, 3))
                }

                function Ss(t, e) {
                    return (yf(t) ? f : kc)(t, Ri(e, 3))
                }

                function ks(t, e, n, r) {
                    t = ou(t) ? t : ma(t), n = n && !r ? Bu(n) : 0;
                    var i = t.length;
                    return n < 0 && (n = Ql(i + n, 0)), Tu(t) ? n <= i && t.indexOf(e, n) > -1 : !!i && S(t, e, n) > -1
                }

                function js(t, e) {
                    return (yf(t) ? v : cr)(t, Ri(e, 3))
                }

                function Os(t, e, n, r) {
                    return null == t ? [] : (yf(e) || (e = null == e ? [] : [e]), n = r ? Q : n, yf(n) || (n = null == n ? [] : [n]), vr(t, e, n))
                }

                function Es(t, e, n) {
                    var r = yf(t) ? y : O,
                        i = arguments.length < 3;
                    return r(t, Ri(e, 4), n, i, Sc)
                }

                function Ms(t, e, n) {
                    var r = yf(t) ? _ : O,
                        i = arguments.length < 3;
                    return r(t, Ri(e, 4), n, i, kc)
                }

                function zs(t, e) {
                    var n = yf(t) ? d : Bn;
                    return e = Ri(e, 3), n(t, function(t, n, r) {
                        return !e(t, n, r)
                    })
                }

                function As(t) {
                    var e = ou(t) ? t : ma(t),
                        n = e.length;
                    return n > 0 ? e[Sr(0, n - 1)] : Q
                }

                function Cs(t, e, n) {
                    var r = -1,
                        i = Fu(t),
                        o = i.length,
                        s = o - 1;
                    for (e = (n ? Yi(t, e, n) : e === Q) ? 1 : _n(Bu(e), 0, o); ++r < e;) {
                        var u = Sr(r, s),
                            a = i[u];
                        i[u] = i[r], i[r] = a
                    }
                    return i.length = e, i
                }

                function Ts(t) {
                    return Cs(t, Mt)
                }

                function Is(t) {
                    if (null == t) return 0;
                    if (ou(t)) {
                        var e = t.length;
                        return e && Tu(t) ? X(t) : e
                    }
                    if (_u(t)) {
                        var n = Di(t);
                        if (n == Ft || n == Ht) return t.size
                    }
                    return ia(t).length
                }

                function Ns(t, e, n) {
                    var r = yf(t) ? b : Mr;
                    return n && Yi(t, e, n) && (e = Q), r(t, Ri(e, 3))
                }

                function Ls() {
                    return _l.now()
                }

                function Rs(t, e) {
                    if ("function" != typeof e) throw new Sl(nt);
                    return t = Bu(t),
                        function() {
                            if (--t < 1) return e.apply(this, arguments)
                        }
                }

                function qs(t, e, n) {
                    return e = n ? Q : e, e = t && null == e ? t.length : e, Mi(t, ht, Q, Q, Q, Q, e)
                }

                function Fs(t, e) {
                    var n;
                    if ("function" != typeof e) throw new Sl(nt);
                    return t = Bu(t),
                        function() {
                            return --t > 0 && (n = e.apply(this, arguments)), t <= 1 && (e = Q), n
                        }
                }

                function Ps(t, e, n) {
                    e = n ? Q : e;
                    var r = Mi(t, at, Q, Q, Q, Q, Q, e);
                    return r.placeholder = Ps.placeholder, r
                }

                function Bs(t, e, n) {
                    e = n ? Q : e;
                    var r = Mi(t, lt, Q, Q, Q, Q, Q, e);
                    return r.placeholder = Bs.placeholder, r
                }

                function Ws(t, e, n) {
                    function r(e) {
                        var n = h,
                            r = d;
                        return h = d = Q, y = e, g = t.apply(r, n)
                    }

                    function i(t) {
                        return y = t, v = $l(u, e), _ ? r(t) : g
                    }

                    function o(t) {
                        var n = t - m,
                            r = t - y,
                            i = e - n;
                        return b ? tc(i, p - r) : i
                    }

                    function s(t) {
                        var n = t - m,
                            r = t - y;
                        return m === Q || n >= e || n < 0 || b && r >= p
                    }

                    function u() {
                        var t = Ls();
                        return s(t) ? a(t) : void(v = $l(u, o(t)))
                    }

                    function a(t) {
                        return v = Q, x && h ? r(t) : (h = d = Q, g)
                    }

                    function l() {
                        y = 0, h = m = d = v = Q
                    }

                    function c() {
                        return v === Q ? g : a(Ls())
                    }

                    function f() {
                        var t = Ls(),
                            n = s(t);
                        if (h = arguments, d = this, m = t, n) {
                            if (v === Q) return i(m);
                            if (b) return v = $l(u, e), r(m)
                        }
                        return v === Q && (v = $l(u, e)), g
                    }
                    var h, d, p, g, v, m, y = 0,
                        _ = !1,
                        b = !1,
                        x = !0;
                    if ("function" != typeof t) throw new Sl(nt);
                    return e = Du(e) || 0, yu(n) && (_ = !!n.leading, b = "maxWait" in n, p = b ? Ql(Du(n.maxWait) || 0, e) : p, x = "trailing" in n ? !!n.trailing : x), f.cancel = l, f.flush = c, f
                }

                function Ds(t) {
                    return Mi(t, pt)
                }

                function Hs(t, e) {
                    if ("function" != typeof t || e && "function" != typeof e) throw new Sl(nt);
                    var n = function r() {
                        var n = arguments,
                            i = e ? e.apply(this, n) : n[0],
                            o = r.cache;
                        if (o.has(i)) return o.get(i);
                        var s = t.apply(this, n);
                        return r.cache = o.set(i, s), s
                    };
                    return n.cache = new(Hs.Cache || Ze), n
                }

                function Us(t) {
                    if ("function" != typeof t) throw new Sl(nt);
                    return function() {
                        return !t.apply(this, arguments)
                    }
                }

                function Vs(t) {
                    return Fs(2, t)
                }

                function $s(t, e) {
                    if ("function" != typeof t) throw new Sl(nt);
                    return e = Ql(e === Q ? t.length - 1 : Bu(e), 0),
                        function() {
                            for (var n = arguments, r = -1, i = Ql(n.length - e, 0), o = Array(i); ++r < i;) o[r] = n[e + r];
                            switch (e) {
                                case 0:
                                    return t.call(this, o);
                                case 1:
                                    return t.call(this, n[0], o);
                                case 2:
                                    return t.call(this, n[0], n[1], o)
                            }
                            var s = Array(e + 1);
                            for (r = -1; ++r < e;) s[r] = n[r];
                            return s[e] = o, a(t, this, s)
                        }
                }

                function Gs(t, e) {
                    if ("function" != typeof t) throw new Sl(nt);
                    return e = e === Q ? 0 : Ql(Bu(e), 0), $s(function(n) {
                        var r = n[e],
                            i = Ur(n, 0, e);
                        return r && m(i, r), a(t, this, i)
                    })
                }

                function Ks(t, e, n) {
                    var r = !0,
                        i = !0;
                    if ("function" != typeof t) throw new Sl(nt);
                    return yu(n) && (r = "leading" in n ? !!n.leading : r, i = "trailing" in n ? !!n.trailing : i), Ws(t, e, {
                        leading: r,
                        maxWait: e,
                        trailing: i
                    })
                }

                function Xs(t) {
                    return qs(t, 1)
                }

                function Zs(t, e) {
                    return e = null == e ? Ga : e, df(e, t)
                }

                function Js() {
                    if (!arguments.length) return [];
                    var t = arguments[0];
                    return yf(t) ? t : [t]
                }

                function Ys(t) {
                    return wn(t, !1, !0)
                }

                function Qs(t, e) {
                    return wn(t, !1, !0, e)
                }

                function tu(t) {
                    return wn(t, !0, !0)
                }

                function eu(t, e) {
                    return wn(t, !0, !0, e)
                }

                function nu(t, e) {
                    return t === e || t !== t && e !== e
                }

                function ru(t) {
                    return su(t) && Al.call(t, "callee") && (!Ul.call(t, "callee") || Il.call(t) == Ct)
                }

                function iu(t) {
                    return _u(t) && Il.call(t) == Kt
                }

                function ou(t) {
                    return null != t && mu(Ac(t)) && !gu(t)
                }

                function su(t) {
                    return _u(t) && ou(t)
                }

                function uu(t) {
                    return t === !0 || t === !1 || _u(t) && Il.call(t) == It
                }

                function au(t) {
                    return _u(t) && Il.call(t) == Nt
                }

                function lu(t) {
                    return !!t && 1 === t.nodeType && _u(t) && !Mu(t)
                }

                function cu(t) {
                    if (ou(t) && (yf(t) || Tu(t) || gu(t.splice) || ru(t) || _f(t))) return !t.length;
                    if (_u(t)) {
                        var e = Di(t);
                        if (e == Ft || e == Ht) return !t.size
                    }
                    for (var n in t)
                        if (Al.call(t, n)) return !1;
                    return !(dc && ia(t).length)
                }

                function fu(t, e) {
                    return nr(t, e)
                }

                function hu(t, e, n) {
                    n = "function" == typeof n ? n : Q;
                    var r = n ? n(t, e) : Q;
                    return r === Q ? nr(t, e, n) : !!r
                }

                function du(t) {
                    return !!_u(t) && (Il.call(t) == Lt || "string" == typeof t.message && "string" == typeof t.name)
                }

                function pu(t) {
                    return "number" == typeof t && Zl(t)
                }

                function gu(t) {
                    var e = yu(t) ? Il.call(t) : "";
                    return e == Rt || e == qt
                }

                function vu(t) {
                    return "number" == typeof t && t == Bu(t)
                }

                function mu(t) {
                    return "number" == typeof t && t > -1 && t % 1 == 0 && t <= jt
                }

                function yu(t) {
                    var e = "undefined" == typeof t ? "undefined" : o(t);
                    return !!t && ("object" == e || "function" == e)
                }

                function _u(t) {
                    return !!t && "object" == ("undefined" == typeof t ? "undefined" : o(t))
                }

                function bu(t) {
                    return _u(t) && Di(t) == Ft
                }

                function xu(t, e) {
                    return t === e || ir(t, e, Fi(e))
                }

                function wu(t, e, n) {
                    return n = "function" == typeof n ? n : Q, ir(t, e, Fi(e), n)
                }

                function Su(t) {
                    return Eu(t) && t != +t
                }

                function ku(t) {
                    if (Tc(t)) throw new bl("This method is not supported with `core-js`. Try https://github.com/es-shims.");
                    return or(t)
                }

                function ju(t) {
                    return null === t
                }

                function Ou(t) {
                    return null == t
                }

                function Eu(t) {
                    return "number" == typeof t || _u(t) && Il.call(t) == Pt
                }

                function Mu(t) {
                    if (!_u(t) || Il.call(t) != Bt || H(t)) return !1;
                    var e = Bi(t);
                    if (null === e) return !0;
                    var n = Al.call(e, "constructor") && e.constructor;
                    return "function" == typeof n && n instanceof n && zl.call(n) == Tl
                }

                function zu(t) {
                    return yu(t) && Il.call(t) == Dt
                }

                function Au(t) {
                    return vu(t) && t >= -jt && t <= jt
                }

                function Cu(t) {
                    return _u(t) && Di(t) == Ht
                }

                function Tu(t) {
                    return "string" == typeof t || !yf(t) && _u(t) && Il.call(t) == Ut
                }

                function Iu(t) {
                    return "symbol" == ("undefined" == typeof t ? "undefined" : o(t)) || _u(t) && Il.call(t) == Vt
                }

                function Nu(t) {
                    return _u(t) && mu(t.length) && !!Mn[Il.call(t)]
                }

                function Lu(t) {
                    return t === Q
                }

                function Ru(t) {
                    return _u(t) && Di(t) == $t
                }

                function qu(t) {
                    return _u(t) && Il.call(t) == Gt
                }

                function Fu(t) {
                    if (!t) return [];
                    if (ou(t)) return Tu(t) ? Z(t) : ri(t);
                    if (Dl && t[Dl]) return U(t[Dl]());
                    var e = Di(t);
                    return (e == Ft ? V : e == Ht ? G : ma)(t)
                }

                function Pu(t) {
                    return t ? (t = Du(t), t === kt || t === -kt ? (t < 0 ? -1 : 1) * Ot : t === t ? t : 0) : 0 === t ? t : 0
                }

                function Bu(t) {
                    var e = Pu(t),
                        n = e % 1;
                    return e === e ? n ? e - n : e : 0
                }

                function Wu(t) {
                    return t ? _n(Bu(t), 0, Mt) : 0
                }

                function Du(t) {
                    if ("number" == typeof t) return t;
                    if (Iu(t)) return Et;
                    if (yu(t)) {
                        var e = gu(t.valueOf) ? t.valueOf() : t;
                        t = yu(e) ? e + "" : e
                    }
                    if ("string" != typeof t) return 0 === t ? t : +t;
                    t = t.replace(be, "");
                    var n = ze.test(t);
                    return n || Ce.test(t) ? Ln(t.slice(2), n ? 2 : 8) : Me.test(t) ? Et : +t
                }

                function Hu(t) {
                    return ii(t, oa(t))
                }

                function Uu(t) {
                    return _n(Bu(t), -jt, jt)
                }

                function Vu(t) {
                    return null == t ? "" : Ir(t)
                }

                function $u(t, e) {
                    var n = Cn(t);
                    return e ? mn(n, e) : n
                }

                function Gu(t, e) {
                    return x(t, Ri(e, 3), Un)
                }

                function Ku(t, e) {
                    return x(t, Ri(e, 3), Vn)
                }

                function Xu(t, e) {
                    return null == t ? t : jc(t, Ri(e, 3), oa)
                }

                function Zu(t, e) {
                    return null == t ? t : Oc(t, Ri(e, 3), oa)
                }

                function Ju(t, e) {
                    return t && Un(t, Ri(e, 3))
                }

                function Yu(t, e) {
                    return t && Vn(t, Ri(e, 3))
                }

                function Qu(t) {
                    return null == t ? [] : $n(t, ia(t))
                }

                function ta(t) {
                    return null == t ? [] : $n(t, oa(t))
                }

                function ea(t, e, n) {
                    var r = null == t ? Q : Gn(t, e);
                    return r === Q ? n : r
                }

                function na(t, e) {
                    return null != t && Ui(t, e, Zn)
                }

                function ra(t, e) {
                    return null != t && Ui(t, e, Jn)
                }

                function ia(t) {
                    var e = ro(t);
                    if (!e && !ou(t)) return ur(t);
                    var n = Ki(t),
                        r = !!n,
                        i = n || [],
                        o = i.length;
                    for (var s in t) !Zn(t, s) || r && ("length" == s || Ji(s, o)) || e && "constructor" == s || i.push(s);
                    return i
                }

                function oa(t) {
                    for (var e = -1, n = ro(t), r = ar(t), i = r.length, o = Ki(t), s = !!o, u = o || [], a = u.length; ++e < i;) {
                        var l = r[e];
                        s && ("length" == l || Ji(l, a)) || "constructor" == l && (n || !Al.call(t, l)) || u.push(l)
                    }
                    return u
                }

                function sa(t, e) {
                    var n = {};
                    return e = Ri(e, 3), Un(t, function(t, r, i) {
                        n[e(t, r, i)] = t
                    }), n
                }

                function ua(t, e) {
                    var n = {};
                    return e = Ri(e, 3), Un(t, function(t, r, i) {
                        n[r] = e(t, r, i)
                    }), n
                }

                function aa(t, e) {
                    return e = Ri(e), yr(t, function(t, n) {
                        return !e(t, n)
                    })
                }

                function la(t, e) {
                    return null == t ? {} : yr(t, Ri(e))
                }

                function ca(t, e, n) {
                    e = Qi(e, t) ? [e] : Hr(e);
                    var r = -1,
                        i = e.length;
                    for (i || (t = Q, i = 1); ++r < i;) {
                        var o = null == t ? Q : t[co(e[r])];
                        o === Q && (r = i, o = n), t = gu(o) ? o.call(t) : o
                    }
                    return t
                }

                function fa(t, e, n) {
                    return null == t ? t : Or(t, e, n)
                }

                function ha(t, e, n, r) {
                    return r = "function" == typeof r ? r : Q, null == t ? t : Or(t, e, n, r)
                }

                function da(t, e, n) {
                    var r = yf(t) || Nu(t);
                    if (e = Ri(e, 4), null == n)
                        if (r || yu(t)) {
                            var i = t.constructor;
                            n = r ? yf(t) ? new i : [] : gu(i) ? Cn(Bi(t)) : {}
                        } else n = {};
                    return (r ? c : Un)(t, function(t, r, i) {
                        return e(n, t, r, i)
                    }), n
                }

                function pa(t, e) {
                    return null == t || Lr(t, e)
                }

                function ga(t, e, n) {
                    return null == t ? t : Rr(t, e, Dr(n))
                }

                function va(t, e, n, r) {
                    return r = "function" == typeof r ? r : Q, null == t ? t : Rr(t, e, Dr(n), r)
                }

                function ma(t) {
                    return t ? T(t, ia(t)) : []
                }

                function ya(t) {
                    return null == t ? [] : T(t, oa(t))
                }

                function _a(t, e, n) {
                    return n === Q && (n = e, e = Q), n !== Q && (n = Du(n), n = n === n ? n : 0), e !== Q && (e = Du(e), e = e === e ? e : 0), _n(Du(t), e, n)
                }

                function ba(t, e, n) {
                    return e = Du(e) || 0, n === Q ? (n = e, e = 0) : n = Du(n) || 0, t = Du(t), Yn(t, e, n)
                }

                function xa(t, e, n) {
                    if (n && "boolean" != typeof n && Yi(t, e, n) && (e = n = Q), n === Q && ("boolean" == typeof e ? (n = e, e = Q) : "boolean" == typeof t && (n = t, t = Q)), t === Q && e === Q ? (t = 0, e = 1) : (t = Du(t) || 0, e === Q ? (e = t, t = 0) : e = Du(e) || 0), t > e) {
                        var r = t;
                        t = e, e = r
                    }
                    if (n || t % 1 || e % 1) {
                        var i = nc();
                        return tc(t + i * (e - t + Nn("1e-" + ((i + "").length - 1))), e)
                    }
                    return Sr(t, e)
                }

                function wa(t) {
                    return Vf(Vu(t).toLowerCase())
                }

                function Sa(t) {
                    return t = Vu(t), t && t.replace(Ie, F).replace(xn, "")
                }

                function ka(t, e, n) {
                    t = Vu(t), e = Ir(e);
                    var r = t.length;
                    return n = n === Q ? r : _n(Bu(n), 0, r), n -= e.length, n >= 0 && t.indexOf(e, n) == n
                }

                function ja(t) {
                    return t = Vu(t), t && fe.test(t) ? t.replace(le, P) : t
                }

                function Oa(t) {
                    return t = Vu(t), t && _e.test(t) ? t.replace(ye, "\\$&") : t
                }

                function Ea(t, e, n) {
                    t = Vu(t), e = Bu(e);
                    var r = e ? X(t) : 0;
                    if (!e || r >= e) return t;
                    var i = (e - r) / 2;
                    return xi(Kl(i), n) + t + xi(Gl(i), n)
                }

                function Ma(t, e, n) {
                    t = Vu(t), e = Bu(e);
                    var r = e ? X(t) : 0;
                    return e && r < e ? t + xi(e - r, n) : t
                }

                function za(t, e, n) {
                    t = Vu(t), e = Bu(e);
                    var r = e ? X(t) : 0;
                    return e && r < e ? xi(e - r, n) + t : t
                }

                function Aa(t, e, n) {
                    return n || null == e ? e = 0 : e && (e = +e), t = Vu(t).replace(be, ""), ec(t, e || (Ee.test(t) ? 16 : 10))
                }

                function Ca(t, e, n) {
                    return e = (n ? Yi(t, e, n) : e === Q) ? 1 : Bu(e), jr(Vu(t), e)
                }

                function Ta() {
                    var t = arguments,
                        e = Vu(t[0]);
                    return t.length < 3 ? e : rc.call(e, t[1], t[2])
                }

                function Ia(t, e, n) {
                    return n && "number" != typeof n && Yi(t, e, n) && (e = n = Q), (n = n === Q ? Mt : n >>> 0) ? (t = Vu(t), t && ("string" == typeof e || null != e && !zu(e)) && (e = Ir(e), "" == e && kn.test(t)) ? Ur(Z(t), 0, n) : oc.call(t, e, n)) : []
                }

                function Na(t, e, n) {
                    return t = Vu(t), n = _n(Bu(n), 0, t.length), t.lastIndexOf(Ir(e), n) == n
                }

                function La(t, n, r) {
                    var i = e.templateSettings;
                    r && Yi(t, n, r) && (n = Q), t = Vu(t), n = kf({}, n, i, hn);
                    var o, s, u = kf({}, n.imports, i.imports, hn),
                        a = ia(u),
                        l = T(u, a),
                        c = 0,
                        f = n.interpolate || Ne,
                        h = "__p += '",
                        d = wl((n.escape || Ne).source + "|" + f.source + "|" + (f === pe ? je : Ne).source + "|" + (n.evaluate || Ne).source + "|$", "g"),
                        p = "//# sourceURL=" + ("sourceURL" in n ? n.sourceURL : "lodash.templateSources[" + ++En + "]") + "\n";
                    t.replace(d, function(e, n, r, i, u, a) {
                        return r || (r = i), h += t.slice(c, a).replace(Le, B), n && (o = !0, h += "' +\n__e(" + n + ") +\n'"), u && (s = !0, h += "';\n" + u + ";\n__p += '"), r && (h += "' +\n((__t = (" + r + ")) == null ? '' : __t) +\n'"), c = a + e.length, e
                    }), h += "';\n";
                    var g = n.variable;
                    g || (h = "with (obj) {\n" + h + "\n}\n"), h = (s ? h.replace(oe, "") : h).replace(se, "$1").replace(ue, "$1;"), h = "function(" + (g || "obj") + ") {\n" + (g ? "" : "obj || (obj = {});\n") + "var __t, __p = ''" + (o ? ", __e = _.escape" : "") + (s ? ", __j = Array.prototype.join;\nfunction print() { __p += __j.call(arguments, '') }\n" : ";\n") + h + "return __p\n}";
                    var v = $f(function() {
                        return Function(a, p + "return " + h).apply(Q, l)
                    });
                    if (v.source = h, du(v)) throw v;
                    return v
                }

                function Ra(t) {
                    return Vu(t).toLowerCase()
                }

                function qa(t) {
                    return Vu(t).toUpperCase()
                }

                function Fa(t, e, n) {
                    if (t = Vu(t), t && (n || e === Q)) return t.replace(be, "");
                    if (!t || !(e = Ir(e))) return t;
                    var r = Z(t),
                        i = Z(e);
                    return Ur(r, N(r, i), L(r, i) + 1).join("")
                }

                function Pa(t, e, n) {
                    if (t = Vu(t), t && (n || e === Q)) return t.replace(we, "");
                    if (!t || !(e = Ir(e))) return t;
                    var r = Z(t);
                    return Ur(r, 0, L(r, Z(e)) + 1).join("")
                }

                function Ba(t, e, n) {
                    if (t = Vu(t), t && (n || e === Q)) return t.replace(xe, "");
                    if (!t || !(e = Ir(e))) return t;
                    var r = Z(t);
                    return Ur(r, N(r, Z(e))).join("")
                }

                function Wa(t, e) {
                    var n = mt,
                        r = yt;
                    if (yu(e)) {
                        var i = "separator" in e ? e.separator : i;
                        n = "length" in e ? Bu(e.length) : n, r = "omission" in e ? Ir(e.omission) : r
                    }
                    t = Vu(t);
                    var o = t.length;
                    if (kn.test(t)) {
                        var s = Z(t);
                        o = s.length
                    }
                    if (n >= o) return t;
                    var u = n - X(r);
                    if (u < 1) return r;
                    var a = s ? Ur(s, 0, u).join("") : t.slice(0, u);
                    if (i === Q) return a + r;
                    if (s && (u += a.length - u), zu(i)) {
                        if (t.slice(u).search(i)) {
                            var l, c = a;
                            for (i.global || (i = wl(i.source, Vu(Oe.exec(i)) + "g")), i.lastIndex = 0; l = i.exec(c);) var f = l.index;
                            a = a.slice(0, f === Q ? u : f)
                        }
                    } else if (t.indexOf(Ir(i), u) != u) {
                        var h = a.lastIndexOf(i);
                        h > -1 && (a = a.slice(0, h))
                    }
                    return a + r
                }

                function Da(t) {
                    return t = Vu(t), t && ce.test(t) ? t.replace(ae, J) : t
                }

                function Ha(t, e, n) {
                    return t = Vu(t), e = n ? Q : e, e === Q && (e = jn.test(t) ? Sn : Se), t.match(e) || []
                }

                function Ua(t) {
                    var e = t ? t.length : 0,
                        n = Ri();
                    return t = e ? v(t, function(t) {
                        if ("function" != typeof t[1]) throw new Sl(nt);
                        return [n(t[0]), t[1]]
                    }) : [], $s(function(n) {
                        for (var r = -1; ++r < e;) {
                            var i = t[r];
                            if (a(i[0], this, n)) return a(i[1], this, n)
                        }
                    })
                }

                function Va(t) {
                    return An(wn(t, !0))
                }

                function $a(t) {
                    return function() {
                        return t
                    }
                }

                function Ga(t) {
                    return t
                }

                function Ka(t) {
                    return sr("function" == typeof t ? t : wn(t, !0))
                }

                function Xa(t) {
                    return fr(wn(t, !0))
                }

                function Za(t, e) {
                    return hr(t, wn(e, !0))
                }

                function Ja(t, e, n) {
                    var r = ia(e),
                        i = $n(e, r);
                    null != n || yu(e) && (i.length || !r.length) || (n = e, e = t, t = this, i = $n(e, ia(e)));
                    var o = !(yu(n) && "chain" in n && !n.chain),
                        s = gu(t);
                    return c(i, function(n) {
                        var r = e[n];
                        t[n] = r, s && (t.prototype[n] = function() {
                            var e = this.__chain__;
                            if (o || e) {
                                var n = t(this.__wrapped__);
                                return (n.__actions__ = ri(this.__actions__)).push({
                                    func: r,
                                    args: arguments,
                                    thisArg: t
                                }), n.__chain__ = e, n
                            }
                            return r.apply(t, m([this.value()], arguments))
                        })
                    }), t
                }

                function Ya() {
                    return Dn._ === this && (Dn._ = Nl), this
                }

                function Qa() {}

                function tl(t) {
                    return t = Bu(t), $s(function(e) {
                        return gr(e, t)
                    })
                }

                function el(t) {
                    return Qi(t) ? _r(co(t)) : br(t)
                }

                function nl(t) {
                    return function(e) {
                        return null == t ? Q : Gn(t, e)
                    }
                }

                function rl() {
                    return []
                }

                function il() {
                    return !1
                }

                function ol() {
                    return {}
                }

                function sl() {
                    return ""
                }

                function ul() {
                    return !0
                }

                function al(t, e) {
                    if (t = Bu(t), t < 1 || t > jt) return [];
                    var n = Mt,
                        r = tc(t, Mt);
                    e = Ri(e), t -= Mt;
                    for (var i = z(r, e); ++n < t;) e(n);
                    return i
                }

                function ll(t) {
                    return yf(t) ? v(t, co) : Iu(t) ? [t] : ri(Nc(t))
                }

                function cl(t) {
                    var e = ++Cl;
                    return Vu(t) + e
                }

                function fl(t) {
                    return t && t.length ? qn(t, Ga, Xn) : Q
                }

                function hl(t, e) {
                    return t && t.length ? qn(t, Ri(e), Xn) : Q
                }

                function dl(t) {
                    return j(t, Ga)
                }

                function pl(t, e) {
                    return j(t, Ri(e))
                }

                function gl(t) {
                    return t && t.length ? qn(t, Ga, lr) : Q
                }

                function vl(t, e) {
                    return t && t.length ? qn(t, Ri(e), lr) : Q
                }

                function ml(t) {
                    return t && t.length ? M(t, Ga) : 0
                }

                function yl(t, e) {
                    return t && t.length ? M(t, Ri(e)) : 0
                }
                t = t ? Hn.defaults({}, t, Hn.pick(Dn, On)) : Dn;
                var _l = t.Date,
                    bl = t.Error,
                    xl = t.Math,
                    wl = t.RegExp,
                    Sl = t.TypeError,
                    kl = t.Array.prototype,
                    jl = t.Object.prototype,
                    Ol = t.String.prototype,
                    El = t["__core-js_shared__"],
                    Ml = function() {
                        var t = /[^.]+$/.exec(El && El.keys && El.keys.IE_PROTO || "");
                        return t ? "Symbol(src)_1." + t : ""
                    }(),
                    zl = t.Function.prototype.toString,
                    Al = jl.hasOwnProperty,
                    Cl = 0,
                    Tl = zl.call(Object),
                    Il = jl.toString,
                    Nl = Dn._,
                    Ll = wl("^" + zl.call(Al).replace(ye, "\\$&").replace(/hasOwnProperty|(function).*?(?=\\\()| for .+?(?=\\\])/g, "$1.*?") + "$"),
                    Rl = Fn ? t.Buffer : Q,
                    ql = t.Reflect,
                    Fl = t.Symbol,
                    Pl = t.Uint8Array,
                    Bl = ql ? ql.enumerate : Q,
                    Wl = Object.getOwnPropertySymbols,
                    Dl = "symbol" == o(Dl = Fl && Fl.iterator) ? Dl : Q,
                    Hl = Object.create,
                    Ul = jl.propertyIsEnumerable,
                    Vl = kl.splice,
                    $l = function(e, n) {
                        return t.setTimeout.call(Dn, e, n)
                    },
                    Gl = xl.ceil,
                    Kl = xl.floor,
                    Xl = Object.getPrototypeOf,
                    Zl = t.isFinite,
                    Jl = kl.join,
                    Yl = Object.keys,
                    Ql = xl.max,
                    tc = xl.min,
                    ec = t.parseInt,
                    nc = xl.random,
                    rc = Ol.replace,
                    ic = kl.reverse,
                    oc = Ol.split,
                    sc = Pi(t, "DataView"),
                    uc = Pi(t, "Map"),
                    ac = Pi(t, "Promise"),
                    lc = Pi(t, "Set"),
                    cc = Pi(t, "WeakMap"),
                    fc = Pi(Object, "create"),
                    hc = cc && new cc,
                    dc = !Ul.call({
                        valueOf: 1
                    }, "valueOf"),
                    pc = {},
                    gc = fo(sc),
                    vc = fo(uc),
                    mc = fo(ac),
                    yc = fo(lc),
                    _c = fo(cc),
                    bc = Fl ? Fl.prototype : Q,
                    xc = bc ? bc.valueOf : Q,
                    wc = bc ? bc.toString : Q;
                e.templateSettings = {
                    escape: he,
                    evaluate: de,
                    interpolate: pe,
                    variable: "",
                    imports: {
                        _: e
                    }
                }, e.prototype = n.prototype, e.prototype.constructor = e, r.prototype = Cn(n.prototype), r.prototype.constructor = r, i.prototype = Cn(n.prototype), i.prototype.constructor = i, Fe.prototype.clear = Pe, Fe.prototype["delete"] = Be, Fe.prototype.get = We, Fe.prototype.has = De, Fe.prototype.set = He, Ue.prototype.clear = Ve, Ue.prototype["delete"] = $e, Ue.prototype.get = Ge, Ue.prototype.has = Ke, Ue.prototype.set = Xe, Ze.prototype.clear = Je, Ze.prototype["delete"] = Ye, Ze.prototype.get = Qe, Ze.prototype.has = tn, Ze.prototype.set = en, nn.prototype.add = nn.prototype.push = rn, nn.prototype.has = on, sn.prototype.clear = un, sn.prototype["delete"] = an, sn.prototype.get = ln, sn.prototype.has = cn, sn.prototype.set = fn;
                var Sc = ai(Un),
                    kc = ai(Vn, !0),
                    jc = li(),
                    Oc = li(!0);
                Bl && !Ul.call({
                    valueOf: 1
                }, "valueOf") && (ar = function(t) {
                    return U(Bl(t))
                });
                var Ec = hc ? function(t, e) {
                        return hc.set(t, e), t
                    } : Ga,
                    Mc = lc && 1 / G(new lc([, -0]))[1] == kt ? function(t) {
                        return new lc(t)
                    } : Qa,
                    zc = hc ? function(t) {
                        return hc.get(t)
                    } : Qa,
                    Ac = _r("length");
                Wl || (Wi = rl);
                var Cc = Wl ? function(t) {
                    for (var e = []; t;) m(e, Wi(t)), t = Bi(t);
                    return e
                } : Wi;
                (sc && Di(new sc(new ArrayBuffer(1))) != Xt || uc && Di(new uc) != Ft || ac && Di(ac.resolve()) != Wt || lc && Di(new lc) != Ht || cc && Di(new cc) != $t) && (Di = function(t) {
                    var e = Il.call(t),
                        n = e == Bt ? t.constructor : Q,
                        r = n ? fo(n) : Q;
                    if (r) switch (r) {
                        case gc:
                            return Xt;
                        case vc:
                            return Ft;
                        case mc:
                            return Wt;
                        case yc:
                            return Ht;
                        case _c:
                            return $t
                    }
                    return e
                });
                var Tc = El ? gu : il,
                    Ic = function() {
                        var t = 0,
                            e = 0;
                        return function(n, r) {
                            var i = Ls(),
                                o = bt - (i - e);
                            if (e = i, o > 0) {
                                if (++t >= _t) return n
                            } else t = 0;
                            return Ec(n, r)
                        }
                    }(),
                    Nc = Hs(function(t) {
                        var e = [];
                        return Vu(t).replace(me, function(t, n, r, i) {
                            e.push(r ? i.replace(ke, "$1") : n || t)
                        }), e
                    }),
                    Lc = $s(function(t, e) {
                        return su(t) ? In(t, Wn(e, 1, su, !0)) : []
                    }),
                    Rc = $s(function(t, e) {
                        var n = To(e);
                        return su(n) && (n = Q), su(t) ? In(t, Wn(e, 1, su, !0), Ri(n)) : []
                    }),
                    qc = $s(function(t, e) {
                        var n = To(e);
                        return su(n) && (n = Q), su(t) ? In(t, Wn(e, 1, su, !0), Q, n) : []
                    }),
                    Fc = $s(function(t) {
                        var e = v(t, Wr);
                        return e.length && e[0] === t[0] ? Qn(e) : []
                    }),
                    Pc = $s(function(t) {
                        var e = To(t),
                            n = v(t, Wr);
                        return e === To(n) ? e = Q : n.pop(), n.length && n[0] === t[0] ? Qn(n, Ri(e)) : []
                    }),
                    Bc = $s(function(t) {
                        var e = To(t),
                            n = v(t, Wr);
                        return e === To(n) ? e = Q : n.pop(), n.length && n[0] === t[0] ? Qn(n, Q, e) : []
                    }),
                    Wc = $s(Lo),
                    Dc = $s(function(t, e) {
                        e = Wn(e, 1);
                        var n = t ? t.length : 0,
                            r = yn(t, e);
                        return wr(t, v(e, function(t) {
                            return Ji(t, n) ? +t : t
                        }).sort(Qr)), r
                    }),
                    Hc = $s(function(t) {
                        return Nr(Wn(t, 1, su, !0))
                    }),
                    Uc = $s(function(t) {
                        var e = To(t);
                        return su(e) && (e = Q), Nr(Wn(t, 1, su, !0), Ri(e))
                    }),
                    Vc = $s(function(t) {
                        var e = To(t);
                        return su(e) && (e = Q), Nr(Wn(t, 1, su, !0), Q, e)
                    }),
                    $c = $s(function(t, e) {
                        return su(t) ? In(t, e) : []
                    }),
                    Gc = $s(function(t) {
                        return Pr(d(t, su))
                    }),
                    Kc = $s(function(t) {
                        var e = To(t);
                        return su(e) && (e = Q), Pr(d(t, su), Ri(e))
                    }),
                    Xc = $s(function(t) {
                        var e = To(t);
                        return su(e) && (e = Q), Pr(d(t, su), Q, e)
                    }),
                    Zc = $s(rs),
                    Jc = $s(function(t) {
                        var e = t.length,
                            n = e > 1 ? t[e - 1] : Q;
                        return n = "function" == typeof n ? (t.pop(), n) : Q, is(t, n)
                    }),
                    Yc = $s(function(t) {
                        t = Wn(t, 1);
                        var e = t.length,
                            n = e ? t[0] : 0,
                            o = this.__wrapped__,
                            s = function(e) {
                                return yn(e, t)
                            };
                        return !(e > 1 || this.__actions__.length) && o instanceof i && Ji(n) ? (o = o.slice(n, +n + (e ? 1 : 0)), o.__actions__.push({
                            func: ls,
                            args: [s],
                            thisArg: Q
                        }), new r(o, this.__chain__).thru(function(t) {
                            return e && !t.length && t.push(Q), t
                        })) : this.thru(s)
                    }),
                    Qc = si(function(t, e, n) {
                        Al.call(t, n) ? ++t[n] : t[n] = 1
                    }),
                    tf = gi(wo),
                    ef = gi(So),
                    nf = si(function(t, e, n) {
                        Al.call(t, n) ? t[n].push(e) : t[n] = [e]
                    }),
                    rf = $s(function(t, e, n) {
                        var r = -1,
                            i = "function" == typeof e,
                            o = Qi(e),
                            s = ou(t) ? Array(t.length) : [];
                        return Sc(t, function(t) {
                            var u = i ? e : o && null != t ? t[e] : Q;
                            s[++r] = u ? a(u, t, n) : er(t, e, n)
                        }), s
                    }),
                    of = si(function(t, e, n) {
                        t[n] = e
                    }),
                    sf = si(function(t, e, n) {
                        t[n ? 0 : 1].push(e)
                    }, function() {
                        return [
                            [],
                            []
                        ]
                    }),
                    uf = $s(function(t, e) {
                        if (null == t) return [];
                        var n = e.length;
                        return n > 1 && Yi(t, e[0], e[1]) ? e = [] : n > 2 && Yi(e[0], e[1], e[2]) && (e = [e[0]]), e = 1 == e.length && yf(e[0]) ? e[0] : Wn(e, 1, Zi), vr(t, e, [])
                    }),
                    af = $s(function(t, e, n) {
                        var r = ot;
                        if (n.length) {
                            var i = $(n, Li(af));
                            r |= ct
                        }
                        return Mi(t, r, e, n, i)
                    }),
                    lf = $s(function(t, e, n) {
                        var r = ot | st;
                        if (n.length) {
                            var i = $(n, Li(lf));
                            r |= ct
                        }
                        return Mi(e, r, t, n, i)
                    }),
                    cf = $s(function(t, e) {
                        return Tn(t, 1, e)
                    }),
                    ff = $s(function(t, e, n) {
                        return Tn(t, Du(e) || 0, n)
                    });
                Hs.Cache = Ze;
                var hf = $s(function(t, e) {
                        e = 1 == e.length && yf(e[0]) ? v(e[0], C(Ri())) : v(Wn(e, 1, Zi), C(Ri()));
                        var n = e.length;
                        return $s(function(r) {
                            for (var i = -1, o = tc(r.length, n); ++i < o;) r[i] = e[i].call(this, r[i]);
                            return a(t, this, r)
                        })
                    }),
                    df = $s(function(t, e) {
                        return Mi(t, ct, Q, e, $(e, Li(df)))
                    }),
                    pf = $s(function(t, e) {
                        return Mi(t, ft, Q, e, $(e, Li(pf)))
                    }),
                    gf = $s(function(t, e) {
                        return Mi(t, dt, Q, Q, Q, Wn(e, 1))
                    }),
                    vf = ki(Xn),
                    mf = ki(function(t, e) {
                        return t >= e
                    }),
                    yf = Array.isArray,
                    _f = Rl ? function(t) {
                        return t instanceof Rl
                    } : il,
                    bf = ki(lr),
                    xf = ki(function(t, e) {
                        return t <= e
                    }),
                    wf = ui(function(t, e) {
                        if (dc || ro(e) || ou(e)) return void ii(e, ia(e), t);
                        for (var n in e) Al.call(e, n) && pn(t, n, e[n])
                    }),
                    Sf = ui(function(t, e) {
                        if (dc || ro(e) || ou(e)) return void ii(e, oa(e), t);
                        for (var n in e) pn(t, n, e[n])
                    }),
                    kf = ui(function(t, e, n, r) {
                        ii(e, oa(e), t, r)
                    }),
                    jf = ui(function(t, e, n, r) {
                        ii(e, ia(e), t, r)
                    }),
                    Of = $s(function(t, e) {
                        return yn(t, Wn(e, 1))
                    }),
                    Ef = $s(function(t) {
                        return t.push(Q, hn), a(kf, Q, t)
                    }),
                    Mf = $s(function(t) {
                        return t.push(Q, uo), a(If, Q, t)
                    }),
                    zf = yi(function(t, e, n) {
                        t[e] = n
                    }, $a(Ga)),
                    Af = yi(function(t, e, n) {
                        Al.call(t, e) ? t[e].push(n) : t[e] = [n]
                    }, Ri),
                    Cf = $s(er),
                    Tf = ui(function(t, e, n) {
                        dr(t, e, n)
                    }),
                    If = ui(function(t, e, n, r) {
                        dr(t, e, n, r)
                    }),
                    Nf = $s(function(t, e) {
                        return null == t ? {} : (e = v(Wn(e, 1), co), mr(t, In(Ii(t), e)))
                    }),
                    Lf = $s(function(t, e) {
                        return null == t ? {} : mr(t, v(Wn(e, 1), co))
                    }),
                    Rf = Ei(ia),
                    qf = Ei(oa),
                    Ff = hi(function(t, e, n) {
                        return e = e.toLowerCase(), t + (n ? wa(e) : e)
                    }),
                    Pf = hi(function(t, e, n) {
                        return t + (n ? "-" : "") + e.toLowerCase()
                    }),
                    Bf = hi(function(t, e, n) {
                        return t + (n ? " " : "") + e.toLowerCase()
                    }),
                    Wf = fi("toLowerCase"),
                    Df = hi(function(t, e, n) {
                        return t + (n ? "_" : "") + e.toLowerCase()
                    }),
                    Hf = hi(function(t, e, n) {
                        return t + (n ? " " : "") + Vf(e)
                    }),
                    Uf = hi(function(t, e, n) {
                        return t + (n ? " " : "") + e.toUpperCase()
                    }),
                    Vf = fi("toUpperCase"),
                    $f = $s(function(t, e) {
                        try {
                            return a(t, Q, e)
                        } catch (n) {
                            return du(n) ? n : new bl(n)
                        }
                    }),
                    Gf = $s(function(t, e) {
                        return c(Wn(e, 1), function(e) {
                            e = co(e), t[e] = af(t[e], t)
                        }), t
                    }),
                    Kf = vi(),
                    Xf = vi(!0),
                    Zf = $s(function(t, e) {
                        return function(n) {
                            return er(n, t, e)
                        }
                    }),
                    Jf = $s(function(t, e) {
                        return function(n) {
                            return er(t, n, e)
                        }
                    }),
                    Yf = bi(v),
                    Qf = bi(h),
                    th = bi(b),
                    eh = Si(),
                    nh = Si(!0),
                    rh = _i(function(t, e) {
                        return t + e
                    }),
                    ih = Oi("ceil"),
                    oh = _i(function(t, e) {
                        return t / e
                    }),
                    sh = Oi("floor"),
                    uh = _i(function(t, e) {
                        return t * e
                    }),
                    ah = Oi("round"),
                    lh = _i(function(t, e) {
                        return t - e
                    });
                return e.after = Rs, e.ary = qs, e.assign = wf, e.assignIn = Sf, e.assignInWith = kf, e.assignWith = jf, e.at = Of, e.before = Fs, e.bind = af, e.bindAll = Gf, e.bindKey = lf, e.castArray = Js, e.chain = us, e.chunk = po, e.compact = go, e.concat = vo, e.cond = Ua, e.conforms = Va, e.constant = $a, e.countBy = Qc, e.create = $u, e.curry = Ps, e.curryRight = Bs, e.debounce = Ws, e.defaults = Ef, e.defaultsDeep = Mf, e.defer = cf, e.delay = ff, e.difference = Lc, e.differenceBy = Rc, e.differenceWith = qc, e.drop = mo, e.dropRight = yo, e.dropRightWhile = _o, e.dropWhile = bo, e.fill = xo, e.filter = ys, e.flatMap = _s, e.flatMapDeep = bs, e.flatMapDepth = xs, e.flatten = ko, e.flattenDeep = jo, e.flattenDepth = Oo, e.flip = Ds, e.flow = Kf, e.flowRight = Xf, e.fromPairs = Eo, e.functions = Qu, e.functionsIn = ta, e.groupBy = nf, e.initial = Ao, e.intersection = Fc, e.intersectionBy = Pc, e.intersectionWith = Bc, e.invert = zf, e.invertBy = Af, e.invokeMap = rf, e.iteratee = Ka, e.keyBy = of , e.keys = ia, e.keysIn = oa, e.map = js, e.mapKeys = sa, e.mapValues = ua, e.matches = Xa, e.matchesProperty = Za, e.memoize = Hs, e.merge = Tf, e.mergeWith = If, e.method = Zf, e.methodOf = Jf, e.mixin = Ja, e.negate = Us, e.nthArg = tl, e.omit = Nf, e.omitBy = aa, e.once = Vs, e.orderBy = Os, e.over = Yf, e.overArgs = hf, e.overEvery = Qf, e.overSome = th, e.partial = df, e.partialRight = pf, e.partition = sf, e.pick = Lf, e.pickBy = la, e.property = el, e.propertyOf = nl, e.pull = Wc, e.pullAll = Lo, e.pullAllBy = Ro, e.pullAllWith = qo, e.pullAt = Dc, e.range = eh, e.rangeRight = nh, e.rearg = gf, e.reject = zs, e.remove = Fo, e.rest = $s, e.reverse = Po, e.sampleSize = Cs, e.set = fa, e.setWith = ha, e.shuffle = Ts, e.slice = Bo, e.sortBy = uf, e.sortedUniq = Go, e.sortedUniqBy = Ko, e.split = Ia, e.spread = Gs, e.tail = Xo, e.take = Zo, e.takeRight = Jo, e.takeRightWhile = Yo, e.takeWhile = Qo, e.tap = as, e.throttle = Ks, e.thru = ls, e.toArray = Fu, e.toPairs = Rf, e.toPairsIn = qf, e.toPath = ll, e.toPlainObject = Hu, e.transform = da, e.unary = Xs, e.union = Hc, e.unionBy = Uc, e.unionWith = Vc, e.uniq = ts, e.uniqBy = es, e.uniqWith = ns, e.unset = pa, e.unzip = rs, e.unzipWith = is, e.update = ga, e.updateWith = va, e.values = ma, e.valuesIn = ya, e.without = $c, e.words = Ha, e.wrap = Zs, e.xor = Gc, e.xorBy = Kc, e.xorWith = Xc, e.zip = Zc, e.zipObject = os, e.zipObjectDeep = ss, e.zipWith = Jc, e.entries = Rf, e.entriesIn = qf, e.extend = Sf, e.extendWith = kf, Ja(e, e), e.add = rh, e.attempt = $f, e.camelCase = Ff, e.capitalize = wa, e.ceil = ih, e.clamp = _a, e.clone = Ys, e.cloneDeep = tu, e.cloneDeepWith = eu, e.cloneWith = Qs, e.deburr = Sa, e.divide = oh, e.endsWith = ka, e.eq = nu, e.escape = ja, e.escapeRegExp = Oa, e.every = ms, e.find = tf, e.findIndex = wo, e.findKey = Gu, e.findLast = ef, e.findLastIndex = So, e.findLastKey = Ku, e.floor = sh, e.forEach = ws, e.forEachRight = Ss, e.forIn = Xu, e.forInRight = Zu, e.forOwn = Ju, e.forOwnRight = Yu, e.get = ea, e.gt = vf, e.gte = mf, e.has = na, e.hasIn = ra, e.head = Mo, e.identity = Ga, e.includes = ks, e.indexOf = zo, e.inRange = ba, e.invoke = Cf, e.isArguments = ru, e.isArray = yf, e.isArrayBuffer = iu, e.isArrayLike = ou, e.isArrayLikeObject = su, e.isBoolean = uu, e.isBuffer = _f, e.isDate = au, e.isElement = lu, e.isEmpty = cu, e.isEqual = fu, e.isEqualWith = hu, e.isError = du, e.isFinite = pu, e.isFunction = gu, e.isInteger = vu, e.isLength = mu, e.isMap = bu, e.isMatch = xu, e.isMatchWith = wu, e.isNaN = Su, e.isNative = ku, e.isNil = Ou, e.isNull = ju, e.isNumber = Eu, e.isObject = yu, e.isObjectLike = _u, e.isPlainObject = Mu, e.isRegExp = zu, e.isSafeInteger = Au, e.isSet = Cu, e.isString = Tu, e.isSymbol = Iu, e.isTypedArray = Nu, e.isUndefined = Lu, e.isWeakMap = Ru, e.isWeakSet = qu, e.join = Co, e.kebabCase = Pf, e.last = To, e.lastIndexOf = Io, e.lowerCase = Bf, e.lowerFirst = Wf, e.lt = bf, e.lte = xf, e.max = fl, e.maxBy = hl, e.mean = dl, e.meanBy = pl, e.min = gl, e.minBy = vl, e.stubArray = rl, e.stubFalse = il, e.stubObject = ol, e.stubString = sl, e.stubTrue = ul, e.multiply = uh, e.nth = No, e.noConflict = Ya, e.noop = Qa, e.now = Ls, e.pad = Ea, e.padEnd = Ma, e.padStart = za, e.parseInt = Aa, e.random = xa, e.reduce = Es, e.reduceRight = Ms, e.repeat = Ca, e.replace = Ta, e.result = ca, e.round = ah, e.runInContext = Y, e.sample = As, e.size = Is, e.snakeCase = Df, e.some = Ns, e.sortedIndex = Wo, e.sortedIndexBy = Do, e.sortedIndexOf = Ho, e.sortedLastIndex = Uo, e.sortedLastIndexBy = Vo, e.sortedLastIndexOf = $o, e.startCase = Hf, e.startsWith = Na, e.subtract = lh, e.sum = ml, e.sumBy = yl, e.template = La, e.times = al, e.toFinite = Pu, e.toInteger = Bu, e.toLength = Wu, e.toLower = Ra, e.toNumber = Du, e.toSafeInteger = Uu, e.toString = Vu, e.toUpper = qa, e.trim = Fa, e.trimEnd = Pa, e.trimStart = Ba, e.truncate = Wa, e.unescape = Da, e.uniqueId = cl, e.upperCase = Uf, e.upperFirst = Vf, e.each = ws, e.eachRight = Ss, e.first = Mo, Ja(e, function() {
                    var t = {};
                    return Un(e, function(n, r) {
                        Al.call(e.prototype, r) || (t[r] = n)
                    }), t
                }(), {
                    chain: !1
                }), e.VERSION = tt, c(["bind", "bindKey", "curry", "curryRight", "partial", "partialRight"], function(t) {
                    e[t].placeholder = e
                }), c(["drop", "take"], function(t, e) {
                    i.prototype[t] = function(n) {
                        var r = this.__filtered__;
                        if (r && !e) return new i(this);
                        n = n === Q ? 1 : Ql(Bu(n), 0);
                        var o = this.clone();
                        return r ? o.__takeCount__ = tc(n, o.__takeCount__) : o.__views__.push({
                            size: tc(n, Mt),
                            type: t + (o.__dir__ < 0 ? "Right" : "")
                        }), o
                    }, i.prototype[t + "Right"] = function(e) {
                        return this.reverse()[t](e).reverse()
                    }
                }), c(["filter", "map", "takeWhile"], function(t, e) {
                    var n = e + 1,
                        r = n == xt || n == St;
                    i.prototype[t] = function(t) {
                        var e = this.clone();
                        return e.__iteratees__.push({
                            iteratee: Ri(t, 3),
                            type: n
                        }), e.__filtered__ = e.__filtered__ || r, e
                    }
                }), c(["head", "last"], function(t, e) {
                    var n = "take" + (e ? "Right" : "");
                    i.prototype[t] = function() {
                        return this[n](1).value()[0]
                    }
                }), c(["initial", "tail"], function(t, e) {
                    var n = "drop" + (e ? "" : "Right");
                    i.prototype[t] = function() {
                        return this.__filtered__ ? new i(this) : this[n](1)
                    }
                }), i.prototype.compact = function() {
                    return this.filter(Ga)
                }, i.prototype.find = function(t) {
                    return this.filter(t).head()
                }, i.prototype.findLast = function(t) {
                    return this.reverse().find(t)
                }, i.prototype.invokeMap = $s(function(t, e) {
                    return "function" == typeof t ? new i(this) : this.map(function(n) {
                        return er(n, t, e)
                    })
                }), i.prototype.reject = function(t) {
                    return t = Ri(t, 3), this.filter(function(e) {
                        return !t(e)
                    })
                }, i.prototype.slice = function(t, e) {
                    t = Bu(t);
                    var n = this;
                    return n.__filtered__ && (t > 0 || e < 0) ? new i(n) : (t < 0 ? n = n.takeRight(-t) : t && (n = n.drop(t)), e !== Q && (e = Bu(e), n = e < 0 ? n.dropRight(-e) : n.take(e - t)), n)
                }, i.prototype.takeRightWhile = function(t) {
                    return this.reverse().takeWhile(t).reverse()
                }, i.prototype.toArray = function() {
                    return this.take(Mt)
                }, Un(i.prototype, function(t, n) {
                    var o = /^(?:filter|find|map|reject)|While$/.test(n),
                        s = /^(?:head|last)$/.test(n),
                        u = e[s ? "take" + ("last" == n ? "Right" : "") : n],
                        a = s || /^find/.test(n);
                    u && (e.prototype[n] = function() {
                        var n = this.__wrapped__,
                            l = s ? [1] : arguments,
                            c = n instanceof i,
                            f = l[0],
                            h = c || yf(n),
                            d = function(t) {
                                var n = u.apply(e, m([t], l));
                                return s && p ? n[0] : n
                            };
                        h && o && "function" == typeof f && 1 != f.length && (c = h = !1);
                        var p = this.__chain__,
                            g = !!this.__actions__.length,
                            v = a && !p,
                            y = c && !g;
                        if (!a && h) {
                            n = y ? n : new i(this);
                            var _ = t.apply(n, l);
                            return _.__actions__.push({
                                func: ls,
                                args: [d],
                                thisArg: Q
                            }), new r(_, p)
                        }
                        return v && y ? t.apply(this, l) : (_ = this.thru(d), v ? s ? _.value()[0] : _.value() : _)
                    })
                }), c(["pop", "push", "shift", "sort", "splice", "unshift"], function(t) {
                    var n = kl[t],
                        r = /^(?:push|sort|unshift)$/.test(t) ? "tap" : "thru",
                        i = /^(?:pop|shift)$/.test(t);
                    e.prototype[t] = function() {
                        var t = arguments;
                        if (i && !this.__chain__) {
                            var e = this.value();
                            return n.apply(yf(e) ? e : [], t)
                        }
                        return this[r](function(e) {
                            return n.apply(yf(e) ? e : [], t)
                        })
                    }
                }), Un(i.prototype, function(t, n) {
                    var r = e[n];
                    if (r) {
                        var i = r.name + "";
                        (pc[i] || (pc[i] = [])).push({
                            name: n,
                            func: r
                        })
                    }
                }), pc[mi(Q, st).name] = [{
                    name: "wrapper",
                    func: Q
                }], i.prototype.clone = R, i.prototype.reverse = Re, i.prototype.value = qe, e.prototype.at = Yc, e.prototype.chain = cs, e.prototype.commit = fs, e.prototype.next = hs, e.prototype.plant = ps, e.prototype.reverse = gs, e.prototype.toJSON = e.prototype.valueOf = e.prototype.value = vs, Dl && (e.prototype[Dl] = ds), e
            }
            var Q, tt = "4.13.1",
                et = 200,
                nt = "Expected a function",
                rt = "__lodash_hash_undefined__",
                it = "__lodash_placeholder__",
                ot = 1,
                st = 2,
                ut = 4,
                at = 8,
                lt = 16,
                ct = 32,
                ft = 64,
                ht = 128,
                dt = 256,
                pt = 512,
                gt = 1,
                vt = 2,
                mt = 30,
                yt = "...",
                _t = 150,
                bt = 16,
                xt = 1,
                wt = 2,
                St = 3,
                kt = 1 / 0,
                jt = 9007199254740991,
                Ot = 1.7976931348623157e308,
                Et = NaN,
                Mt = 4294967295,
                zt = Mt - 1,
                At = Mt >>> 1,
                Ct = "[object Arguments]",
                Tt = "[object Array]",
                It = "[object Boolean]",
                Nt = "[object Date]",
                Lt = "[object Error]",
                Rt = "[object Function]",
                qt = "[object GeneratorFunction]",
                Ft = "[object Map]",
                Pt = "[object Number]",
                Bt = "[object Object]",
                Wt = "[object Promise]",
                Dt = "[object RegExp]",
                Ht = "[object Set]",
                Ut = "[object String]",
                Vt = "[object Symbol]",
                $t = "[object WeakMap]",
                Gt = "[object WeakSet]",
                Kt = "[object ArrayBuffer]",
                Xt = "[object DataView]",
                Zt = "[object Float32Array]",
                Jt = "[object Float64Array]",
                Yt = "[object Int8Array]",
                Qt = "[object Int16Array]",
                te = "[object Int32Array]",
                ee = "[object Uint8Array]",
                ne = "[object Uint8ClampedArray]",
                re = "[object Uint16Array]",
                ie = "[object Uint32Array]",
                oe = /\b__p \+= '';/g,
                se = /\b(__p \+=) '' \+/g,
                ue = /(__e\(.*?\)|\b__t\)) \+\n'';/g,
                ae = /&(?:amp|lt|gt|quot|#39|#96);/g,
                le = /[&<>"'`]/g,
                ce = RegExp(ae.source),
                fe = RegExp(le.source),
                he = /<%-([\s\S]+?)%>/g,
                de = /<%([\s\S]+?)%>/g,
                pe = /<%=([\s\S]+?)%>/g,
                ge = /\.|\[(?:[^[\]]*|(["'])(?:(?!\1)[^\\]|\\.)*?\1)\]/,
                ve = /^\w*$/,
                me = /[^.[\]]+|\[(?:(-?\d+(?:\.\d+)?)|(["'])((?:(?!\2)[^\\]|\\.)*?)\2)\]|(?=(\.|\[\])(?:\4|$))/g,
                ye = /[\\^$.*+?()[\]{}|]/g,
                _e = RegExp(ye.source),
                be = /^\s+|\s+$/g,
                xe = /^\s+/,
                we = /\s+$/,
                Se = /[a-zA-Z0-9]+/g,
                ke = /\\(\\)?/g,
                je = /\$\{([^\\}]*(?:\\.[^\\}]*)*)\}/g,
                Oe = /\w*$/,
                Ee = /^0x/i,
                Me = /^[-+]0x[0-9a-f]+$/i,
                ze = /^0b[01]+$/i,
                Ae = /^\[object .+?Constructor\]$/,
                Ce = /^0o[0-7]+$/i,
                Te = /^(?:0|[1-9]\d*)$/,
                Ie = /[\xc0-\xd6\xd8-\xde\xdf-\xf6\xf8-\xff]/g,
                Ne = /($^)/,
                Le = /['\n\r\u2028\u2029\\]/g,
                Re = "\\ud800-\\udfff",
                qe = "\\u0300-\\u036f\\ufe20-\\ufe23",
                Fe = "\\u20d0-\\u20f0",
                Pe = "\\u2700-\\u27bf",
                Be = "a-z\\xdf-\\xf6\\xf8-\\xff",
                We = "\\xac\\xb1\\xd7\\xf7",
                De = "\\x00-\\x2f\\x3a-\\x40\\x5b-\\x60\\x7b-\\xbf",
                He = "\\u2000-\\u206f",
                Ue = " \\t\\x0b\\f\\xa0\\ufeff\\n\\r\\u2028\\u2029\\u1680\\u180e\\u2000\\u2001\\u2002\\u2003\\u2004\\u2005\\u2006\\u2007\\u2008\\u2009\\u200a\\u202f\\u205f\\u3000",
                Ve = "A-Z\\xc0-\\xd6\\xd8-\\xde",
                $e = "\\ufe0e\\ufe0f",
                Ge = We + De + He + Ue,
                Ke = "['’]",
                Xe = "[" + Re + "]",
                Ze = "[" + Ge + "]",
                Je = "[" + qe + Fe + "]",
                Ye = "\\d+",
                Qe = "[" + Pe + "]",
                tn = "[" + Be + "]",
                en = "[^" + Re + Ge + Ye + Pe + Be + Ve + "]",
                nn = "\\ud83c[\\udffb-\\udfff]",
                rn = "(?:" + Je + "|" + nn + ")",
                on = "[^" + Re + "]",
                sn = "(?:\\ud83c[\\udde6-\\uddff]){2}",
                un = "[\\ud800-\\udbff][\\udc00-\\udfff]",
                an = "[" + Ve + "]",
                ln = "\\u200d",
                cn = "(?:" + tn + "|" + en + ")",
                fn = "(?:" + an + "|" + en + ")",
                hn = "(?:" + Ke + "(?:d|ll|m|re|s|t|ve))?",
                dn = "(?:" + Ke + "(?:D|LL|M|RE|S|T|VE))?",
                pn = rn + "?",
                gn = "[" + $e + "]?",
                vn = "(?:" + ln + "(?:" + [on, sn, un].join("|") + ")" + gn + pn + ")*",
                mn = gn + pn + vn,
                yn = "(?:" + [Qe, sn, un].join("|") + ")" + mn,
                _n = "(?:" + [on + Je + "?", Je, sn, un, Xe].join("|") + ")",
                bn = RegExp(Ke, "g"),
                xn = RegExp(Je, "g"),
                wn = RegExp(nn + "(?=" + nn + ")|" + _n + mn, "g"),
                Sn = RegExp([an + "?" + tn + "+" + hn + "(?=" + [Ze, an, "$"].join("|") + ")", fn + "+" + dn + "(?=" + [Ze, an + cn, "$"].join("|") + ")", an + "?" + cn + "+" + hn, an + "+" + dn, Ye, yn].join("|"), "g"),
                kn = RegExp("[" + ln + Re + qe + Fe + $e + "]"),
                jn = /[a-z][A-Z]|[A-Z]{2,}[a-z]|[0-9][a-zA-Z]|[a-zA-Z][0-9]|[^a-zA-Z0-9 ]/,
                On = ["Array", "Buffer", "DataView", "Date", "Error", "Float32Array", "Float64Array", "Function", "Int8Array", "Int16Array", "Int32Array", "Map", "Math", "Object", "Promise", "Reflect", "RegExp", "Set", "String", "Symbol", "TypeError", "Uint8Array", "Uint8ClampedArray", "Uint16Array", "Uint32Array", "WeakMap", "_", "isFinite", "parseInt", "setTimeout"],
                En = -1,
                Mn = {};
            Mn[Zt] = Mn[Jt] = Mn[Yt] = Mn[Qt] = Mn[te] = Mn[ee] = Mn[ne] = Mn[re] = Mn[ie] = !0, Mn[Ct] = Mn[Tt] = Mn[Kt] = Mn[It] = Mn[Xt] = Mn[Nt] = Mn[Lt] = Mn[Rt] = Mn[Ft] = Mn[Pt] = Mn[Bt] = Mn[Dt] = Mn[Ht] = Mn[Ut] = Mn[$t] = !1;
            var zn = {};
            zn[Ct] = zn[Tt] = zn[Kt] = zn[Xt] = zn[It] = zn[Nt] = zn[Zt] = zn[Jt] = zn[Yt] = zn[Qt] = zn[te] = zn[Ft] = zn[Pt] = zn[Bt] = zn[Dt] = zn[Ht] = zn[Ut] = zn[Vt] = zn[ee] = zn[ne] = zn[re] = zn[ie] = !0, zn[Lt] = zn[Rt] = zn[$t] = !1;
            var An = {
                    "À": "A",
                    "Á": "A",
                    "Â": "A",
                    "Ã": "A",
                    "Ä": "A",
                    "Å": "A",
                    "à": "a",
                    "á": "a",
                    "â": "a",
                    "ã": "a",
                    "ä": "a",
                    "å": "a",
                    "Ç": "C",
                    "ç": "c",
                    "Ð": "D",
                    "ð": "d",
                    "È": "E",
                    "É": "E",
                    "Ê": "E",
                    "Ë": "E",
                    "è": "e",
                    "é": "e",
                    "ê": "e",
                    "ë": "e",
                    "Ì": "I",
                    "Í": "I",
                    "Î": "I",
                    "Ï": "I",
                    "ì": "i",
                    "í": "i",
                    "î": "i",
                    "ï": "i",
                    "Ñ": "N",
                    "ñ": "n",
                    "Ò": "O",
                    "Ó": "O",
                    "Ô": "O",
                    "Õ": "O",
                    "Ö": "O",
                    "Ø": "O",
                    "ò": "o",
                    "ó": "o",
                    "ô": "o",
                    "õ": "o",
                    "ö": "o",
                    "ø": "o",
                    "Ù": "U",
                    "Ú": "U",
                    "Û": "U",
                    "Ü": "U",
                    "ù": "u",
                    "ú": "u",
                    "û": "u",
                    "ü": "u",
                    "Ý": "Y",
                    "ý": "y",
                    "ÿ": "y",
                    "Æ": "Ae",
                    "æ": "ae",
                    "Þ": "Th",
                    "þ": "th",
                    "ß": "ss"
                },
                Cn = {
                    "&": "&amp;",
                    "<": "&lt;",
                    ">": "&gt;",
                    '"': "&quot;",
                    "'": "&#39;",
                    "`": "&#96;"
                },
                Tn = {
                    "&amp;": "&",
                    "&lt;": "<",
                    "&gt;": ">",
                    "&quot;": '"',
                    "&#39;": "'",
                    "&#96;": "`"
                },
                In = {
                    "\\": "\\",
                    "'": "'",
                    "\n": "n",
                    "\r": "r",
                    "\u2028": "u2028",
                    "\u2029": "u2029"
                },
                Nn = parseFloat,
                Ln = parseInt,
                Rn = "object" == o(e) && e,
                qn = Rn && "object" == o(t) && t,
                Fn = qn && qn.exports === Rn,
                Pn = R("object" == ("undefined" == typeof i ? "undefined" : o(i)) && i),
                Bn = R("object" == ("undefined" == typeof self ? "undefined" : o(self)) && self),
                Wn = R("object" == o(this) && this),
                Dn = Pn || Bn || Wn || Function("return this")(),
                Hn = Y();
            (Bn || {})._ = Hn, "object" == o(n(53)) && n(53) ? (r = function() {
                return Hn
            }.call(e, n, e, t), !(r !== Q && (t.exports = r))) : qn ? ((qn.exports = Hn)._ = Hn, Rn._ = Hn) : Dn._ = Hn
        }).call(void 0)
    }).call(e, n(20)(t), function() {
        return this
    }())
}, function(t, e, n) {
    "use strict";

    function r(t) {
        var e = t.parentNode,
            n = t.nextSibling;
        return e.removeChild(t),
            function() {
                n ? e.insertBefore(t, n) : e.appendChild(t)
            }
    }
    var i, o = n(54),
        s = n(2),
        u = function(t) {
            for (var e; e = t.lastChild;) t.removeChild(e)
        };
    t.exports = i = s.extend({
        renderSubviews: function() {
            var t = this.el,
                e = void 0 != t.parentNode;
            if (e) var n = r(t);
            u(t);
            for (var i, s, a = document.createDocumentFragment(), l = this._views(), c = o.sortBy(l, function(t) {
                    return t.ordering
                }), f = 0; f < c.length; f++) i = c[f], i.render(), s = i.el, null != s && a.appendChild(s);
            return t.appendChild(a), e && n(), t
        },
        addView: function(t, e) {
            var n = this._views();
            if (null == e) throw "Invalid plugin. ";
            return null == e.ordering && (e.ordering = t), n[t] = e
        },
        removeViews: function() {
            var t, e, n = this._views();
            for (e in n) t = n[e], t.undelegateEvents(), t.unbind(), null != t.removeViews && t.removeViews(), t.remove();
            return this.views = {}
        },
        removeView: function(t) {
            var e = this._views();
            return e[t].remove(), delete e[t]
        },
        getView: function(t) {
            return this._views()[t]
        },
        remove: function() {
            return this.removeViews(), s.prototype.remove.apply(this)
        },
        _views: function() {
            return null == this.views && (this.views = {}), this.views
        }
    })
}, function(t, e, n) {
    var r;
    (function(t) {
        "use strict";
        var i = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        };
        /*!
         * jBone v1.2.0 - 2016-04-13 - Library for DOM manipulation
         *
         * http://jbone.js.org
         *
         * Copyright 2016 Alexey Kupriyanenko
         * Released under the MIT license.
         */
        ! function(o) {
            function s(t) {
                var e = t.length,
                    n = "undefined" == typeof t ? "undefined" : i(t);
                return !_(n) && t !== o && (!(1 !== t.nodeType || !e) || b(n) || 0 === e || "number" == typeof e && e > 0 && e - 1 in t)
            }

            function u(t, e) {
                var n, r;
                this.originalEvent = t, r = function(t, e) {
                    "preventDefault" === t ? this[t] = function() {
                        return this.defaultPrevented = !0, e[t]()
                    } : "stopImmediatePropagation" === t ? this[t] = function() {
                        return this.immediatePropagationStopped = !0, e[t]()
                    } : _(e[t]) ? this[t] = function() {
                        return e[t]()
                    } : this[t] = e[t]
                };
                for (n in t)(t[n] || "function" == typeof t[n]) && r.call(this, n, t);
                x.extend(this, e, {
                    isImmediatePropagationStopped: function() {
                        return !!this.immediatePropagationStopped
                    }
                })
            }
            var a, l = o.$,
                c = o.jBone,
                f = /^<(\w+)\s*\/?>$/,
                h = /^(?:[^#<]*(<[\w\W]+>)[^>]*$|#([\w\-]*)$)/,
                d = [].slice,
                p = [].splice,
                g = Object.keys,
                v = document,
                m = function(t) {
                    return "string" == typeof t
                },
                y = function(t) {
                    return t instanceof Object
                },
                _ = function(t) {
                    return "[object Function]" === {}.toString.call(t)
                },
                b = function(t) {
                    return Array.isArray(t)
                },
                x = function(t, e) {
                    return new a.init(t, e)
                };
            x.noConflict = function() {
                return o.$ = l, o.jBone = c, x
            }, a = x.fn = x.prototype = {
                init: function(t, e) {
                    var n, r, i, o;
                    if (!t) return this;
                    if (m(t)) {
                        if (r = f.exec(t)) return this[0] = v.createElement(r[1]), this.length = 1, y(e) && this.attr(e), this;
                        if ((r = h.exec(t)) && r[1]) {
                            for (o = v.createDocumentFragment(), i = v.createElement("div"), i.innerHTML = t; i.lastChild;) o.appendChild(i.firstChild);
                            return n = d.call(o.childNodes), x.merge(this, n)
                        }
                        if (x.isElement(e)) return x(e).find(t);
                        try {
                            return n = v.querySelectorAll(t), x.merge(this, n)
                        } catch (s) {
                            return this
                        }
                    }
                    return t.nodeType ? (this[0] = t, this.length = 1, this) : _(t) ? t() : t instanceof x ? t : x.makeArray(t, this)
                },
                pop: [].pop,
                push: [].push,
                reverse: [].reverse,
                shift: [].shift,
                sort: [].sort,
                splice: [].splice,
                slice: [].slice,
                indexOf: [].indexOf,
                forEach: [].forEach,
                unshift: [].unshift,
                concat: [].concat,
                join: [].join,
                every: [].every,
                some: [].some,
                filter: [].filter,
                map: [].map,
                reduce: [].reduce,
                reduceRight: [].reduceRight,
                length: 0
            }, a.constructor = x, a.init.prototype = a, x.setId = function(t) {
                var e = t.jid;
                t === o ? e = "window" : void 0 === t.jid && (t.jid = e = ++x._cache.jid), x._cache.events[e] || (x._cache.events[e] = {})
            }, x.getData = function(t) {
                t = t instanceof x ? t[0] : t;
                var e = t === o ? "window" : t.jid;
                return {
                    jid: e,
                    events: x._cache.events[e]
                }
            }, x.isElement = function(t) {
                return t && t instanceof x || t instanceof HTMLElement || m(t)
            }, x._cache = {
                events: {},
                jid: 0
            }, a.pushStack = function(t) {
                return x.merge(this.constructor(), t)
            }, x.merge = function(t, e) {
                for (var n = e.length, r = t.length, i = 0; i < n;) t[r++] = e[i++];
                return t.length = r, t
            }, x.contains = function(t, e) {
                return t.contains(e)
            }, x.extend = function(t) {
                var e;
                return p.call(arguments, 1).forEach(function(n) {
                    if (e = t, n)
                        for (var r in n) e[r] = n[r]
                }), t
            }, x.makeArray = function(t, e) {
                var n = e || [];
                return null !== t && (s(t) ? x.merge(n, m(t) ? [t] : t) : n.push(t)), n
            }, x.unique = function(t) {
                if (null == t) return [];
                for (var e = [], n = 0, r = t.length; n < r; n++) {
                    var i = t[n];
                    e.indexOf(i) < 0 && e.push(i)
                }
                return e
            }, x.Event = function(t, e) {
                var n, r;
                return t.type && !e && (e = t, t = t.type), n = t.split(".").splice(1).join("."), r = t.split(".")[0], t = v.createEvent("Event"), t.initEvent(r, !0, !0), x.extend(t, {
                    namespace: n,
                    isDefaultPrevented: function() {
                        return t.defaultPrevented
                    }
                }, e)
            }, x.event = {
                add: function(t, e, n, r, i) {
                    x.setId(t);
                    var o, s, u, a = function(e) {
                            x.event.dispatch.call(t, e)
                        },
                        l = x.getData(t).events;
                    for (e = e.split(" "), s = e.length; s--;) u = e[s], o = u.split(".")[0], l[o] = l[o] || [], l[o].length ? a = l[o][0].fn : t.addEventListener && t.addEventListener(o, a, !1), l[o].push({
                        namespace: u.split(".").splice(1).join("."),
                        fn: a,
                        selector: i,
                        data: r,
                        originfn: n
                    })
                },
                remove: function(t, e, n, r) {
                    var i, o, s = function(t, e, r, i, o) {
                            var s;
                            (n && o.originfn === n || !n) && (s = o.fn), t[e][r].fn === s && (t[e].splice(r, 1), t[e].length || i.removeEventListener(e, s))
                        },
                        u = x.getData(t).events;
                    if (u) return !e && u ? g(u).forEach(function(e) {
                        for (o = u[e], i = o.length; i--;) s(u, e, i, t, o[i])
                    }) : void e.split(" ").forEach(function(e) {
                        var n, a = e.split(".")[0],
                            l = e.split(".").splice(1).join(".");
                        if (u[a])
                            for (o = u[a], i = o.length; i--;) n = o[i], (!l || l && n.namespace === l) && (!r || r && n.selector === r) && s(u, a, i, t, n);
                        else l && g(u).forEach(function(e) {
                            for (o = u[e], i = o.length; i--;) n = o[i], n.namespace.split(".")[0] === l.split(".")[0] && s(u, e, i, t, n)
                        })
                    })
                },
                trigger: function(t, e) {
                    var n = [];
                    m(e) ? n = e.split(" ").map(function(t) {
                        return x.Event(t)
                    }) : (e = e instanceof Event ? e : x.Event(e), n = [e]), n.forEach(function(e) {
                        e.type && t.dispatchEvent && t.dispatchEvent(e)
                    })
                },
                dispatch: function(t) {
                    for (var e, n, r, i, o, s = 0, a = 0, l = this, c = x.getData(l).events[t.type], f = c.length, h = [], d = []; s < f; s++) h.push(c[s]);
                    for (s = 0, f = h.length; s < f && ~c.indexOf(h[s]) && (!i || !i.isImmediatePropagationStopped()); s++)
                        if (n = null, o = {}, r = h[s], r.data && (o.data = r.data), r.selector) {
                            if (~(d = x(l).find(r.selector)).indexOf(t.target) && (n = t.target) || l !== t.target && l.contains(t.target)) {
                                if (!n)
                                    for (e = d.length, a = 0; a < e; a++) d[a] && d[a].contains(t.target) && (n = d[a]);
                                if (!n) continue;
                                o.currentTarget = n, i = new u(t, o), t.namespace && t.namespace !== r.namespace || r.originfn.call(n, i)
                            }
                        } else i = new u(t, o), t.namespace && t.namespace !== r.namespace || r.originfn.call(l, i)
                }
            }, a.on = function(t, e, n, r) {
                var i = this.length,
                    o = 0;
                if (null == n && null == r ? (r = e, n = e = void 0) : null == r && ("string" == typeof e ? (r = n, n = void 0) : (r = n, n = e, e = void 0)), !r) return this;
                for (; o < i; o++) x.event.add(this[o], t, r, n, e);
                return this
            }, a.one = function(t) {
                var e, n = arguments,
                    r = 0,
                    i = this.length,
                    o = d.call(n, 1, n.length - 1),
                    s = d.call(n, -1)[0];
                for (e = function(e) {
                        var n = x(e);
                        t.split(" ").forEach(function(t) {
                            var r = function i(r) {
                                n.off(t, i), s.call(e, r)
                            };
                            n.on.apply(n, [t].concat(o, r))
                        })
                    }; r < i; r++) e(this[r]);
                return this
            }, a.trigger = function(t) {
                var e = 0,
                    n = this.length;
                if (!t) return this;
                for (; e < n; e++) x.event.trigger(this[e], t);
                return this
            }, a.off = function(t, e, n) {
                var r = 0,
                    i = this.length;
                for (_(e) && (n = e, e = void 0); r < i; r++) x.event.remove(this[r], t, n, e);
                return this
            }, a.find = function(t) {
                for (var e = [], n = 0, r = this.length, i = function(n) {
                        _(n.querySelectorAll) && [].forEach.call(n.querySelectorAll(t), function(t) {
                            e.push(t)
                        })
                    }; n < r; n++) i(this[n]);
                return x(e)
            }, a.get = function(t) {
                return null != t ? t < 0 ? this[t + this.length] : this[t] : d.call(this)
            }, a.eq = function(t) {
                return x(this[t])
            }, a.parent = function() {
                for (var t, e = [], n = 0, r = this.length; n < r; n++) !~e.indexOf(t = this[n].parentElement) && t && e.push(t);
                return x(e)
            }, a.toArray = function() {
                return d.call(this)
            }, a.is = function() {
                var t = arguments;
                return this.some(function(e) {
                    return e.tagName.toLowerCase() === t[0]
                })
            }, a.has = function() {
                var t = arguments;
                return this.some(function(e) {
                    return e.querySelectorAll(t[0]).length
                })
            }, a.add = function(t, e) {
                return this.pushStack(x.unique(x.merge(this.get(), x(t, e))))
            }, a.attr = function(t, e) {
                var n, r = arguments,
                    i = 0,
                    o = this.length;
                if (m(t) && 1 === r.length) return this[0] && this[0].getAttribute(t);
                for (2 === r.length ? n = function(n) {
                        n.setAttribute(t, e)
                    } : y(t) && (n = function(e) {
                        g(t).forEach(function(n) {
                            e.setAttribute(n, t[n])
                        })
                    }); i < o; i++) n(this[i]);
                return this
            }, a.removeAttr = function(t) {
                for (var e = 0, n = this.length; e < n; e++) this[e].removeAttribute(t);
                return this
            }, a.val = function(t) {
                var e = 0,
                    n = this.length;
                if (0 === arguments.length) return this[0] && this[0].value;
                for (; e < n; e++) this[e].value = t;
                return this
            }, a.css = function(t, e) {
                var n, r = arguments,
                    i = 0,
                    s = this.length;
                if (m(t) && 1 === r.length) return this[0] && o.getComputedStyle(this[0])[t];
                for (2 === r.length ? n = function(n) {
                        n.style[t] = e
                    } : y(t) && (n = function(e) {
                        g(t).forEach(function(n) {
                            e.style[n] = t[n]
                        })
                    }); i < s; i++) n(this[i]);
                return this
            }, a.data = function(t, e) {
                var n, r = arguments,
                    i = {},
                    o = 0,
                    s = this.length,
                    u = function(t, e, n) {
                        y(n) ? (t.jdata = t.jdata || {}, t.jdata[e] = n) : t.dataset[e] = n
                    },
                    a = function(t) {
                        return "true" === t || "false" !== t && t
                    };
                if (0 === r.length) return this[0].jdata && (i = this[0].jdata), g(this[0].dataset).forEach(function(t) {
                    i[t] = a(this[0].dataset[t])
                }, this), i;
                if (1 === r.length && m(t)) return this[0] && a(this[0].dataset[t] || this[0].jdata && this[0].jdata[t]);
                for (1 === r.length && y(t) ? n = function(e) {
                        g(t).forEach(function(n) {
                            u(e, n, t[n])
                        })
                    } : 2 === r.length && (n = function(n) {
                        u(n, t, e)
                    }); o < s; o++) n(this[o]);
                return this
            }, a.removeData = function(t) {
                for (var e, n, r = 0, i = this.length; r < i; r++)
                    if (e = this[r].jdata, n = this[r].dataset, t) e && e[t] && delete e[t], delete n[t];
                    else {
                        for (t in e) delete e[t];
                        for (t in n) delete n[t]
                    } return this
            }, a.addClass = function(t) {
                for (var e = 0, n = 0, r = this.length, i = t ? t.trim().split(/\s+/) : []; e < r; e++)
                    for (n = 0, n = 0; n < i.length; n++) this[e].classList.add(i[n]);
                return this
            }, a.removeClass = function(t) {
                for (var e = 0, n = 0, r = this.length, i = t ? t.trim().split(/\s+/) : []; e < r; e++)
                    for (n = 0, n = 0; n < i.length; n++) this[e].classList.remove(i[n]);
                return this
            }, a.toggleClass = function(t, e) {
                var n = 0,
                    r = this.length,
                    i = "toggle";
                if (e === !0 && (i = "add") || e === !1 && (i = "remove"), t)
                    for (; n < r; n++) this[n].classList[i](t);
                return this
            }, a.hasClass = function(t) {
                var e = 0,
                    n = this.length;
                if (t)
                    for (; e < n; e++)
                        if (this[e].classList.contains(t)) return !0;
                return !1
            }, a.html = function(t) {
                var e, n = arguments;
                return 1 === n.length && void 0 !== t ? this.empty().append(t) : 0 === n.length && (e = this[0]) ? e.innerHTML : this
            }, a.append = function(t) {
                var e, n = 0,
                    r = this.length;
                for (m(t) && h.exec(t) ? t = x(t) : y(t) || (t = document.createTextNode(t)), t = t instanceof x ? t : x(t), e = function(e, n) {
                        t.forEach(function(t) {
                            n ? e.appendChild(t.cloneNode(!0)) : e.appendChild(t)
                        })
                    }; n < r; n++) e(this[n], n);
                return this
            }, a.appendTo = function(t) {
                return x(t).append(this), this
            }, a.empty = function() {
                for (var t, e = 0, n = this.length; e < n; e++)
                    for (t = this[e]; t.lastChild;) t.removeChild(t.lastChild);
                return this
            }, a.remove = function() {
                var t, e = 0,
                    n = this.length;
                for (this.off(); e < n; e++) t = this[e], delete t.jdata, t.parentNode && t.parentNode.removeChild(t);
                return this
            }, "object" === i(t) && t && "object" === i(t.exports) ? t.exports = x : (r = function() {
                return x
            }.call(e, n, e, t), !(void 0 !== r && (t.exports = r)), o.jBone = o.$ = x)
        }(window)
    }).call(e, n(20)(t))
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(91),
        i = r.extend({
            buildDOM: function() {
                return this.on("new:node", this.buildNode), this.on("new:button", this.buildButton), this.on("new:menu", this.buildMenu), r.prototype.buildDOM.call(this)
            },
            buildNode: function(t) {
                if (null != this.g) return t.style.lineHeight = this.g.menuconfig.get("menuItemLineHeight")
            },
            buildButton: function(t) {
                if (null != this.g) return t.style.fontSize = this.g.menuconfig.get("menuFontsize"), t.style.marginLeft = this.g.menuconfig.get("menuMarginLeft"), t.style.padding = this.g.menuconfig.get("menuPadding")
            },
            buildMenu: function(t) {
                if (null != this.g) return t.style.fontSize = this.g.menuconfig.get("menuItemFontsize")
            }
        });
    e["default"] = i
}, function(t, e) {
    "use strict";
    var n = {};
    n.removeToInsertLater = function(t) {
        var e, n;
        return n = t.parentNode, e = t.nextSibling, n.removeChild(t),
            function() {
                e ? n.insertBefore(t, e) : n.appendChild(t)
            }
    }, n.removeAllChilds = function(t) {
        var e;
        for (e = 0; t.firstChild;) e++, t.removeChild(t.firstChild)
    }, t.exports = n
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e.columnsel = e.rowsel = e.possel = e.sel = void 0;
    var r = n(3),
        i = n(1).Model,
        o = i.extend({
            defaults: {
                type: "super"
            }
        }),
        s = o.extend({
            defaults: (0, r.extend)({}, o.prototype.defaults, {
                type: "row",
                seqId: ""
            }),
            inRow: function(t) {
                return t === this.get("seqId")
            },
            inColumn: function(t) {
                return !0
            },
            getLength: function() {
                return 1
            }
        }),
        u = o.extend({
            defaults: (0, r.extend)({}, o.prototype.defaults, {
                type: "column",
                xStart: -1,
                xEnd: -1
            }),
            inRow: function() {
                return !0
            },
            inColumn: function(t) {
                return xStart <= t && t <= xEnd
            },
            getLength: function() {
                return xEnd - xStart
            }
        }),
        a = s.extend((0, r.extend)({}, (0, r.pick)(u, "inColumn"), (0, r.pick)(u, "getLength"), {
            defaults: (0, r.extend)({}, u.prototype.defaults, s.prototype.defaults, {
                type: "pos"
            })
        }));
    e.sel = o, e.possel = a, e.rowsel = s, e.columnsel = u
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(8),
        i = n(3),
        o = n(1).Collection,
        s = o.extend({
            model: r.sel,
            initialize: function(t, e) {
                if ("undefined" != typeof e && null !== e) return this.g = e.g, this.listenTo(this.g, "residue:click", function(t) {
                    return this._handleE(t.evt, new r.possel({
                        xStart: t.rowPos,
                        xEnd: t.rowPos,
                        seqId: t.seqId
                    }))
                }), this.listenTo(this.g, "row:click", function(t) {
                    return this._handleE(t.evt, new r.rowsel({
                        seqId: t.seqId
                    }))
                }), this.listenTo(this.g, "column:click", function(t) {
                    return this._handleE(t.evt, new r.columnsel({
                        xStart: t.rowPos,
                        xEnd: t.rowPos + t.stepSize - 1
                    }))
                })
            },
            getSelForRow: function(t) {
                return this.filter(function(e) {
                    return e.inRow(t)
                })
            },
            getSelForColumns: function(t) {
                return this.filter(function(e) {
                    return e.inColumn(t)
                })
            },
            addJSON: function(t) {
                return this.add(this._fromJSON(t))
            },
            _fromJSON: function(t) {
                switch (t.type) {
                    case "column":
                        return new r.columnsel(t);
                    case "row":
                        return new r.rowsel(t);
                    case "pos":
                        return new r.possel(t)
                }
            },
            resetJSON: function(t) {
                return t = t.map(this._fromJSON), this.reset(t)
            },
            getBlocksForRow: function(t, e) {
                for (var n, r = this.filter(function(e) {
                        return e.inRow(t)
                    }), i = [], o = function(t, n) {
                        var n = r[t];
                        return "row" === n.attributes.type ? (i = function() {
                            var t = [],
                                n = 0;
                            if (0 <= e)
                                for (; n <= e;) t.push(n++);
                            else
                                for (; n >= e;) t.push(n--);
                            return t
                        }(), "break") : void(i = i.concat(function() {
                            var t = [],
                                e = n.attributes.xStart;
                            if (n.attributes.xStart <= n.attributes.xEnd)
                                for (; e <= n.attributes.xEnd;) t.push(e++);
                            else
                                for (; e >= n.attributes.xEnd;) t.push(e--);
                            return t
                        }()))
                    }, s = 0; s < r.length && "break" !== o(s, n); s++);
                return i
            },
            getAllColumnBlocks: function(t) {
                var e = (t.maxLen, t.withPos, []),
                    n = void 0;
                n = t.withPos ? this.filter(function(t) {
                    return null != t.get("xStart")
                }) : this.filter(function(t) {
                    return "column" === t.get("type")
                });
                for (var r, o = function(t, r) {
                        var r = n[t];
                        e = e.concat(function() {
                            var t = [],
                                e = r.attributes.xStart;
                            if (r.attributes.xStart <= r.attributes.xEnd)
                                for (; e <= r.attributes.xEnd;) t.push(e++);
                            else
                                for (; e >= r.attributes.xEnd;) t.push(e--);
                            return t
                        }())
                    }, s = 0; s < n.length; s++) o(s, r);
                return e = (0, i.uniq)(e)
            },
            invertRow: function(t) {
                var e = this.where({
                    type: "row"
                });
                e = e.map(function(t) {
                    return t.attributes.seqId
                });
                for (var n, o = (0, i.filter)(t, function(t) {
                        return !(e.indexOf(t) >= 0)
                    }), s = [], u = 0; u < o.length; u++) {
                    var n = o[u];
                    s.push(new r.rowsel({
                        seqId: n
                    }))
                }
                return this.reset(s)
            },
            invertCol: function(t) {
                var e = this.where({
                        type: "column"
                    }).reduce(function(t, e) {
                        return t.concat(function() {
                            var t = [],
                                n = e.attributes.xStart;
                            if (e.attributes.xStart <= e.attributes.xEnd)
                                for (; n <= e.attributes.xEnd;) t.push(n++);
                            else
                                for (; n >= e.attributes.xEnd;) t.push(n--);
                            return t
                        }())
                    }, []),
                    n = (0, i.filter)(t, function(t) {
                        return !(e.indexOf(t) >= 0)
                    });
                if (0 !== n.length) {
                    for (var o, s = [], u = n[0], a = u, l = 0; l < n.length; l++) o = n[l], a + 1 === o ? a = o : (s.push(new r.columnsel({
                        xStart: u,
                        xEnd: a
                    })), u = a = o);
                    return u !== a && s.push(new r.columnsel({
                        xStart: u,
                        xEnd: n[n.length - 1]
                    })), this.reset(s)
                }
            },
            _handleE: function(t, e) {
                return t.ctrlKey || t.metaKey ? this.add(e) : this.reset([e])
            },
            _reduceColumns: function() {
                return this.each(function(t, e, n) {
                    for (var r = (0, i.filter)(n, function(t) {
                            return "column" === t.get("type")
                        }), o = t.get("xStart"), s = t.get("xEnd"), u = (0, i.filter)(r, function(t) {
                            return t.get("xEnd") === o - 1
                        }), a = 0; a < u.length; a++) u[a].set("xEnd", o);
                    for (var l = (0, i.filter)(r, function(t) {
                            return t.get("xStart") === s + 1
                        }), c = 0; c < l.length; c++) l[c].set("xStart", s);
                    if (u.length > 0 || l.length > 0) return t.collection.remove(t)
                })
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e.version = e.xhr = e.seqs = e.parser = e.newick = e.matrix = e.gff = e.fasta = e.clustal = void 0;
    var i = n(61),
        o = r(i),
        s = n(62),
        u = r(s),
        a = n(63),
        l = r(a),
        c = n(65),
        f = r(c),
        h = n(67),
        d = r(h),
        p = n(11),
        g = r(p),
        v = n(24),
        m = r(v);
    e.clustal = o["default"], e.fasta = u["default"], e.gff = l["default"], e.matrix = f["default"], e.newick = d["default"], e.parser = g["default"], e.seqs = m["default"];
    var y = n(21);
    e.xhr = y;
    var _ = "imported";
    "undefined" != typeof IO_VERSION && (_ = IO_VERSION), e.version = _
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        },
        i = n(21),
        o = {};
    e["default"] = o, o.read = function(t, e) {
        var n = function(t) {
            return function(n, r, i) {
                return o._onRetrieval(n, i, e, t)
            }
        }(this);
        return "undefined" == typeof e ? new Promise(function(r, o) {
            e = function(t, e) {
                t ? o(t) : r(e)
            }, i(t, n)
        }) : i(t, n)
    }, o._onRetrieval = function(t, e, n, r) {
        var i = void 0;
        return "undefined" != typeof t && (i = r.parse(e)), n.call(r, t, i)
    }, o.extend = function(t, e) {
        return extend(o, t, e)
    }, o.mixin = function(t) {
        var e = ["read"];
        return "object" !== ("undefined" == typeof t ? "undefined" : r(t)) && (t = t.prototype), e.forEach(function(e) {
            t[e] = o[e]
        }, this), t
    }
}, function(t, e) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var n = "http://www.w3.org/2000/svg",
        r = function(t, e) {
            for (var n in e) {
                var r = e[n];
                t.setAttributeNS(null, n, r)
            }
            return t
        },
        i = function(t) {
            var e = document.createElementNS(n, "svg");
            return e.setAttribute("width", t.width), e.setAttribute("height", t.height), e
        },
        o = function(t) {
            return r(document.createElementNS(n, "rect"), t)
        },
        s = function(t) {
            return r(document.createElementNS(n, "line"), t)
        },
        u = function(t) {
            return r(document.createElementNS(n, "polygon"), t)
        };
    e.base = i, e.line = s, e.rect = o, e.polygon = u
}, function(t, e, n) {
    "use strict";
    var r = n(71);
    r.onAll = function(t, e) {
        return this.on("all", t, e), this
    }, r.oldMixin = r.mixin, r.mixin = function(t) {
        r.oldMixin(t);
        for (var e = ["onAll"], n = 0; n < e.length; n++) {
            var i = e[n];
            t[i] = this[i]
        }
        return t
    }, t.exports = r
}, function(t, e, n) {
    (function(t) {
        "use strict";
        var e = {};
        e.d = e.defaultValue = function(t, e) {
            return void 0 === t ? "function" == typeof t ? e() : e : t
        }, e.id = function(t) {
            return document.getElementById(t)
        }, e.mk = function(t) {
            return document.createElement(t)
        }, void 0 !== t && void 0 !== t.exports && (t.exports = e)
    }).call(e, n(20)(t))
}, function(t, e) {
    "use strict";

    function n(t, e) {
        if (!(t instanceof e)) throw new TypeError("Cannot call a class as a function")
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = function() {
            function t(t, e) {
                for (var n = 0; n < e.length; n++) {
                    var r = e[n];
                    r.enumerable = r.enumerable || !1, r.configurable = !0, "value" in r && (r.writable = !0), Object.defineProperty(t, r.key, r)
                }
            }
            return function(e, n, r) {
                return n && t(e.prototype, n), r && t(e, r), e
            }
        }(),
        i = function() {
            function t() {
                n(this, t)
            }
            return r(t, null, [{
                key: "randomInt",
                value: function(t, e) {
                    if ("undefined" == typeof e || null === e) {
                        var n = [0, t];
                        t = n[0], e = n[1]
                    }
                    if (t > e) {
                        var r = [e, t];
                        t = r[0], e = r[1]
                    }
                    return Math.floor(Math.random() * (e - t + 1) + t)
                }
            }, {
                key: "uniqueId",
                value: function() {
                    for (var t = arguments.length <= 0 || void 0 === arguments[0] ? 8 : arguments[0], e = ""; e.length < t;) e += Math.random().toString(36).substr(2);
                    return e.substr(0, t)
                }
            }, {
                key: "getRandomInt",
                value: function(t, e) {
                    return Math.floor(Math.random() * (e - t + 1)) + t
                }
            }]), t
        }();
    e["default"] = i
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(10),
        i = n(3),
        o = n(87),
        s = n(88),
        u = {
            openInJalview: function(t, e) {
                "." === t.charAt(0) && (t = document.URL.substr(0, document.URL.lastIndexOf("/")) + "/" + t), t.indexOf("http") < 0 && (t = "http://" + window.location.hostname + t), t = encodeURIComponent(t);
                var n = "http://www.jalview.org/services/launchApp?open=" + t;
                return n += "&colour=" + e, window.open(n, "_blank")
            },
            publishWeb: function(t, e) {
                var n = r.fasta.write(t.seqs.toJSON());
                return n = encodeURIComponent(n), (0, r.xhr)({
                    method: "POST",
                    body: "sprunge=" + n,
                    uri: "http://sprunge.biojs.net",
                    headers: {
                        "Content-Type": "application/x-www-form-urlencoded"
                    }
                }, function(t, n, r) {
                    return e(r.trim())
                })
            },
            shareLink: function(t, e) {
                var n = t.g.config.get("importURL"),
                    r = "http://msa.biojs.net/app/?seq=",
                    i = function(t) {
                        var n = r + t;
                        if (e) return e(n)
                    };
                return n ? i(n) : u.publishWeb(t, i)
            },
            saveAsFile: function(t, e) {
                var n = r.fasta.write(t.seqs.toJSON());
                return s(new Blob([n], {
                    type: "text/plain"
                }), e)
            },
            saveSelection: function(t, e) {
                var n = t.g.selcol.pluck("seqId");
                if (n.length > 0) {
                    n = t.seqs.filter(function(t) {
                        return n.indexOf(t.get("id")) >= 0
                    });
                    for (var i = n.length - 1, o = 0; 0 < i ? o <= i : o >= i; 0 < i ? o++ : o--) n[o] = n[o].toJSON()
                } else n = t.seqs.toJSON(), console.warn("no selection found");
                var u = r.fasta.write(n);
                return s(new Blob([u], {
                    type: "text/plain"
                }), e)
            },
            saveAnnots: function(t, e) {
                var n = t.seqs.map(function(t) {
                    if (n = t.get("features"), 0 !== n.length) {
                        var e = t.get("name");
                        return n.each(function(t) {
                            return t.set("seqname", e)
                        }), n.toJSON()
                    }
                });
                n = (0, i.flatten)((0, i.compact)(n));
                var o = r.gff.exportLines(n);
                return s(new Blob([o], {
                    type: "text/plain"
                }), e)
            },
            saveAsImg: function(t, e) {
                var n = t.getView("stage").getView("body").getView("seqblock").el;
                if ("undefined" != typeof n && null !== n) return s(o(n.toDataURL("image/png")), e, "image/png")
            }
        };
    e["default"] = u
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = function() {
            function t(t, e) {
                var n = [],
                    r = !0,
                    i = !1,
                    o = void 0;
                try {
                    for (var s, u = t[Symbol.iterator](); !(r = (s = u.next()).done) && (n.push(s.value), !e || n.length !== e); r = !0);
                } catch (a) {
                    i = !0, o = a
                } finally {
                    try {
                        !r && u["return"] && u["return"]()
                    } finally {
                        if (i) throw o
                    }
                }
                return n
            }
            return function(e, n) {
                if (Array.isArray(e)) return e;
                if (Symbol.iterator in Object(e)) return t(e, n);
                throw new TypeError("Invalid attempt to destructure non-iterable instance")
            }
        }(),
        o = n(3),
        s = n(10),
        u = n(38),
        a = r(u),
        l = function(t) {
            return this.msa = t, this
        },
        c = {
            guessFileFromText: function(t, e) {
                if ("undefined" == typeof t || null === t) return console.warn("invalid file format"), ["", "error"];
                switch ((0, a["default"])(t, e)) {
                    case "clustal":
                        var n = s.clustal,
                            r = "seqs";
                        break;
                    case "fasta":
                        n = s.fasta, r = "seqs";
                        break;
                    case "newick":
                        r = "newick";
                        break;
                    case "gff":
                        n = s.gff, r = "features";
                        break;
                    default:
                        alert("Unknown file format. Please contact us on Github for help.")
                }
                return [n, r]
            },
            parseText: function(t, e) {
                var n = this.guessFileFromText(t, e),
                    r = i(n, 2),
                    o = r[0],
                    s = r[1];
                return "seqs" === s ? [o.parse(t), s] : "features" === s ? [o.parseSeqs(t), s] : [t, s]
            },
            importFiles: function(t) {
                var e = this;
                return function() {
                    for (var n = [], r = t.length - 1, i = 0; 0 < r ? i <= r : i >= r; 0 < r ? i++ : i--) {
                        var o = t[i],
                            s = new FileReader;
                        s.onload = function(t) {
                            return e.importFile(t.target.result)
                        }, n.push(s.readAsText(o))
                    }
                    return n
                }()
            },
            importFile: function(t, e) {
                var n = this;
                e = e || {}, e.name = t.name;
                var r, o = this.parseText(t, e),
                    s = i(o, 2),
                    u = s[0],
                    a = s[1];
                return "error" === a ? (alert("An error happened"), "error") : ("seqs" === a ? (this.msa.seqs.reset(u), this.msa.g.config.set("url", "userimport"), this.msa.g.trigger("url:userImport")) : "features" === a ? this.msa.seqs.addFeatures(u) : "newick" === a ? this.msa.u.tree.loadTree(function() {
                    return n.msa.u.tree.showTree(t)
                }) : alert("Unknown file!"), r = t.name)
            },
            importURL: function(t, e) {
                var n = this;
                return t = this.msa.u.proxy.corsURL(t), this.msa.g.config.set("url", t), (0, s.xhr)({
                    url: t,
                    timeout: 0
                }, function(r, i, o) {
                    return r ? console.error(r) : "error" !== n.importFile(o, {
                        url: t
                    }) ? (n.msa.g.trigger("import:url", t), e ? e() : void 0) : void 0
                })
            }
        };
    (0, o.extend)(l.prototype, c), e["default"] = l
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(14),
        i = {
            loadScript: function(t, e) {
                var n = r.mk("script");
                return n.type = "text/javascript", n.src = t, n.async = !0, n.onload = n.onreadystatechange = function() {
                    if (!(t || this.readyState && "complete" !== this.readyState)) {
                        var t = !0;
                        return e()
                    }
                }, document.getElementsByTagName("script")[0].parentNode.appendChild(n)
            },
            joinCb: function(t, e, n) {
                e = e || 1;
                var r = 0,
                    i = function(t, e) {
                        return "undefined" == typeof t || null === t ? o() : function() {
                            var n;
                            return n = "apply", t.indexOf(n) >= 0 && t.apply(e, arguments), o()
                        }
                    },
                    o = function() {
                        if (r++, r === e) return t.call(n)
                    };
                return i
            }
        };
    e["default"] = i
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(3),
        i = function(t) {
            return this.g = t.g, this
        },
        o = {
            corsURL: function(t) {
                return document.URL.indexOf("localhost") >= 0 && "/" === t[0] ? t : "." === t.charAt(0) || "/" === t.charAt(0) ? t : (this.g.config.get("importProxyStripHttp") && (t = t.replace("http://", ""), t = t.replace("https://", "")), t = this.g.config.get("importProxy") + t)
            }
        };
    (0, r.extend)(i.prototype, o), e["default"] = i
}, function(t, e) {
    "use strict";
    t.exports = function(t) {
        return t.webpackPolyfill || (t.deprecate = function() {}, t.paths = [], t.children = [], t.webpackPolyfill = 1), t
    }
}, function(t, e, n) {
    "use strict";

    function r(t, e) {
        for (var n = 0; n < t.length; n++) e(t[n])
    }

    function i(t) {
        for (var e in t)
            if (t.hasOwnProperty(e)) return !1;
        return !0
    }

    function o(t, e, n) {
        var r = t;
        return f(e) ? (n = e, "string" == typeof t && (r = {
            uri: t
        })) : r = d(e, {
            uri: t
        }), r.callback = n, r
    }

    function s(t, e, n) {
        return e = o(t, e, n), u(e)
    }

    function u(t) {
        function e() {
            4 === f.readyState && o()
        }

        function n() {
            var t = void 0;
            if (t = f.response ? f.response : f.responseText || a(f), x) try {
                t = JSON.parse(t)
            } catch (e) {}
            return t
        }

        function r(t) {
            clearTimeout(g), t instanceof Error || (t = new Error("" + (t || "Unknown XMLHttpRequest Error"))), t.statusCode = 0, u(t, c), u = l
        }

        function o() {
            if (!p) {
                var e;
                clearTimeout(g), e = t.useXDR && void 0 === f.status ? 200 : 1223 === f.status ? 204 : f.status;
                var r = c,
                    i = null;
                0 !== e ? (r = {
                    body: n(),
                    statusCode: e,
                    method: m,
                    headers: {},
                    url: v,
                    rawRequest: f
                }, f.getAllResponseHeaders && (r.headers = h(f.getAllResponseHeaders()))) : i = new Error("Internal XMLHttpRequest Error"), u(i, r, r.body), u = l
            }
        }
        var u = t.callback;
        if ("undefined" == typeof u) throw new Error("callback argument missing");
        var c = {
                body: void 0,
                headers: {},
                statusCode: 0,
                method: m,
                url: v,
                rawRequest: f
            },
            f = t.xhr || null;
        f || (f = t.cors || t.useXDR ? new s.XDomainRequest : new s.XMLHttpRequest);
        var d, p, g, v = f.url = t.uri || t.url,
            m = f.method = t.method || "GET",
            y = t.body || t.data || null,
            _ = f.headers = t.headers || {},
            b = !!t.sync,
            x = !1;
        if ("json" in t && (x = !0, _.accept || _.Accept || (_.Accept = "application/json"), "GET" !== m && "HEAD" !== m && (_["content-type"] || _["Content-Type"] || (_["Content-Type"] = "application/json"), y = JSON.stringify(t.json))), f.onreadystatechange = e, f.onload = o, f.onerror = r, f.onprogress = function() {}, f.ontimeout = r, f.open(m, v, !b, t.username, t.password), b || (f.withCredentials = !!t.withCredentials), !b && t.timeout > 0 && (g = setTimeout(function() {
                p = !0, f.abort("timeout");
                var t = new Error("XMLHttpRequest timeout");
                t.code = "ETIMEDOUT", r(t)
            }, t.timeout)), f.setRequestHeader)
            for (d in _) _.hasOwnProperty(d) && f.setRequestHeader(d, _[d]);
        else if (t.headers && !i(t.headers)) throw new Error("Headers cannot be set on an XDomainRequest object");
        return "responseType" in t && (f.responseType = t.responseType), "beforeSend" in t && "function" == typeof t.beforeSend && t.beforeSend(f), f.send(y), f
    }

    function a(t) {
        if ("document" === t.responseType) return t.responseXML;
        var e = 204 === t.status && t.responseXML && "parsererror" === t.responseXML.documentElement.nodeName;
        return "" !== t.responseType || e ? null : t.responseXML
    }

    function l() {}
    var c = n(147),
        f = n(52),
        h = n(150),
        d = n(151);
    t.exports = s, s.XMLHttpRequest = c.XMLHttpRequest || l, s.XDomainRequest = "withCredentials" in new s.XMLHttpRequest ? s.XMLHttpRequest : c.XDomainRequest, r(["get", "put", "post", "patch", "head", "delete"], function(t) {
        s["delete" === t ? "del" : t] = function(e, n, r) {
            return n = o(e, n, r), n.method = t.toUpperCase(), u(n)
        }
    })
}, function(t, e, n) {
    "use strict";
    t.exports = n(56)
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    ! function(s) {
        "object" === o(e) ? t.exports = s() : (r = s, i = "function" == typeof r ? r.call(e, n, e, t) : r, !(void 0 !== i && (t.exports = i)))
    }(function() {
        "use strict";
        var t = {
            has: function(t, e) {
                return Object.prototype.hasOwnProperty.call(t, e)
            },
            extend: function(t) {
                for (var e = 1; e < arguments.length; ++e) {
                    var n = arguments[e];
                    if (n)
                        for (var r in n) t[r] = n[r]
                }
                return t
            }
        };
        return function(e, n) {
            var r, i = this;
            r = e && t.has(e, "constructor") ? e.constructor : function() {
                return i.apply(this, arguments)
            }, t.extend(r, i, n);
            var o = function() {
                this.constructor = r
            };
            return o.prototype = i.prototype, r.prototype = new o, e && t.extend(r.prototype, e), r.__super__ = i.prototype, r
        }
    })
}, function(t, e) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var n = {};
    e["default"] = n, n.getMeta = function(t) {
        var e, n, r = !1,
            i = !1,
            o = {},
            s = {},
            u = t.split(" ");
        if (u.length >= 1 ? (r = u.shift(), i = u.join(" ")) : r = t, r) {
            var a = r.split("|");
            for (e = a.pop(), s.en = e; 0 != a.length;) {
                var l = a.shift(),
                    c = a.shift();
                o[l] = c
            }
        } else e = r;
        if (i) {
            var f = i.split("=");
            if (f.length > 1) {
                var h, d;
                f.length - 1, f.forEach(function(t) {
                    t = t.trim();
                    var e, r = t.split(" ");
                    r.length > 1 ? (d = r.pop(), e = r.join(" ")) : e = t, h ? s[h.toLowerCase()] = e : n = e, h = d
                })
            } else n = f.shift()
        }
        var p = {
            name: e,
            ids: o,
            details: s
        };
        return n && (p.desc = n), p
    };
    var r = {
        sp: {
            link: "http://www.uniprot.org/%s",
            name: "Uniprot"
        },
        tr: {
            link: "http://www.uniprot.org/%s",
            name: "Trembl"
        },
        gb: {
            link: "http://www.ncbi.nlm.nih.gov/nuccore/%s",
            name: "Genbank"
        },
        pdb: {
            link: "http://www.rcsb.org/pdb/explore/explore.do?structureId=%s",
            name: "PDB"
        }
    };
    n.buildLinks = function(t) {
        var e = {};
        return t = t || {}, Object.keys(t).forEach(function(n) {
            if (n in r) {
                var i = r[n],
                    o = i.link.replace("%s", t[n]);
                e[i.name] = o
            }
        }), e
    }, n.contains = function(t, e) {
        return "".indexOf.call(t, e, 0) !== -1
    }, n.splitNChars = function(t, e) {
        var n, r;
        e = e || 80;
        var i = [];
        for (n = 0, r = t.length - 1; n <= r; n += e) i.push(t.substr(n, e));
        return i
    }, n.reverse = function(t) {
        return t.split("").reverse().join("")
    }, n.complement = function(t) {
        var e = t + "",
            n = [
                [/g/g, "0"],
                [/c/g, "1"],
                [/0/g, "c"],
                [/1/g, "g"],
                [/G/g, "0"],
                [/C/g, "1"],
                [/0/g, "C"],
                [/1/g, "G"],
                [/a/g, "0"],
                [/t/g, "1"],
                [/0/g, "t"],
                [/1/g, "a"],
                [/A/g, "0"],
                [/T/g, "1"],
                [/0/g, "T"],
                [/1/g, "A"]
            ];
        for (var r in n) e = e.replace(n[r][0], n[r][1]);
        return e
    }, n.reverseComplement = function(t) {
        return n.reverse(n.complement(t))
    }, n.model = function(t, e, n) {
        this.seq = t, this.name = e, this.id = n, this.ids = {}
    }
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    n(90), t.exports = r = i.extend({
        constructor: function(t, e) {
            return this.g = e.g, i.apply(this, arguments), this
        },
        defaults: {
            currentSize: 5,
            step: 1,
            originalSize: !1,
            scaleCategories: [{
                columnWidth: 1,
                markerStepSize: 20,
                stepSize: 0
            }, {
                columnWidth: 3,
                markerStepSize: 20,
                stepSize: 0
            }, {
                columnWidth: 5,
                markerStepSize: 10,
                stepSize: 0
            }, {
                columnWidth: 9,
                markerStepSize: 5,
                stepSize: 1
            }, {
                columnWidth: 15,
                markerStepSize: 2,
                stepSize: 1
            }, {
                columnWidth: 20,
                markerStepSize: 1,
                stepSize: 1
            }, {
                columnWidth: 30,
                markerStepSize: 1,
                stepSize: 1
            }, {
                columnWidth: 45,
                markerStepSize: 1,
                stepSize: 1
            }]
        },
        initialize: function(t) {
            var e = this.get("scaleCategories"),
                n = this.g.zoomer.get("columnWidth") || this._getScaleInfo().columnWidth,
                r = _.find(e, {
                    columnWidth: n
                });
            if (!r) {
                var i = this._insertScaleCategory(n);
                r = e[i], this.set("currentSize", i + 1)
            }
            var o = this.get("currentSize");
            return this.set("originalSize", o), this.setSize(o), this
        },
        _insertScaleCategory: function(t) {
            var e = this.get("scaleCategories"),
                n = _.findLastIndex(e, function(e) {
                    return e.columnWidth < t
                }),
                r = e[n],
                i = n + 1,
                o = {
                    columnWidth: t,
                    markerStepSize: r.markerStepSize,
                    stepSize: r.markerStepSize
                };
            return e.splice(i, 0, o), this.set("scaleCategories", e), i
        },
        getSizeRange: function() {
            return [1, this.get("scaleCategories").length]
        },
        bigger: function() {
            return this.setSize(this.get("currentSize") + this.get("step"))
        },
        smaller: function() {
            return this.setSize(this.get("currentSize") - this.get("step"))
        },
        reset: function() {
            return this.setSize(this.get("originalSize"))
        },
        setSize: function(t) {
            var e = this.getSizeRange();
            t = parseInt(t), t = t < e[0] ? e[0] : t > e[1] ? e[1] : t, this.set("currentSize", t);
            var n = this._getScaleInfo();
            return this.g.zoomer.set({
                columnWidth: n.columnWidth,
                stepSize: n.stepSize,
                markerStepSize: n.markerStepSize
            }), this
        },
        getSize: function() {
            return this.get("currentSize")
        },
        _getScaleInfo: function() {
            var t = this.getSize(),
                e = this.get("scaleCategories");
            return t > 0 && t <= e.length ? e[t - 1] : void console.error("out of bounds")
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(99),
        o = n(1).Model;
    t.exports = r = o.extend({
        defaults: {
            scheme: "taylor",
            colorBackground: !0,
            showLowerCase: !0,
            opacity: .6
        },
        initialize: function(t, e, n) {
            return this.colors = new i({
                seqs: e,
                conservation: function() {
                    return n.scale(n.conservation())
                }
            }), n.on("reset", function() {
                if ("dyn" === this.getSelectedScheme().type) {
                    var t;
                    if (t = "reset", this.getSelectedScheme().indexOf(t) >= 0) return this.getSelectedScheme().reset()
                }
            }, this)
        },
        addStaticScheme: function(t, e) {
            return this.colors.addStaticScheme(t, e)
        },
        addDynScheme: function(t, e) {
            return this.colors.addDynScheme(t, e)
        },
        getScheme: function(t) {
            return this.colors.getScheme(t)
        },
        getSelectedScheme: function() {
            return this.colors.getScheme(this.get("scheme"))
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    t.exports = r = i.extend({
        initialize: function(t, e) {
            return null == this.get("hidden") && this.set("hidden", []), this.stats = e
        },
        calcHiddenColumns: function(t) {
            for (var e, n = this.get("hidden"), r = t, i = 0; i < n.length; i++) e = n[i], e <= r && r++;
            return r - t
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    t.exports = r = i.extend({
        defaults: {
            registerMouseHover: !1,
            registerMouseClicks: !0,
            importProxy: "https://cors-anywhere.herokuapp.com/",
            importProxyStripHttp: !0,
            eventBus: !0,
            alphabetSize: 20,
            dropImport: !1,
            debug: !1,
            hasRef: !1,
            bootstrapMenu: !1,
            manualRendering: !1
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(18),
        o = n(1).Model;
    t.exports = r = o.extend({
        initialize: function(t) {
            return this.g = t
        },
        development: {
            "msa-tnt": "/node_modules/msa-tnt/build/bundle.js",
            "biojs-io-newick": "/node_modules/biojs-io-newick/build/biojs-io-newick.min.js"
        },
        loadPackage: function(t, e) {
            try {
                return e(n(155)(t))
            } catch (r) {
                return i["default"].loadScript(this._pkgURL(t), e)
            }
        },
        loadPackages: function(t, e) {
            var n = this,
                r = i["default"].joinCb(function() {
                    return e()
                }, t.length);
            return t.forEach(function(t) {
                return n.loadPackage(t, r)
            })
        },
        _pkgURL: function(t) {
            if (this.g.config.get("debug")) var e = this.development[t];
            else e = "http://wzrd.in/bundle/" + t + "@latest";
            return e
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    t.exports = r = i.extend({
        defaults: {
            searchText: ""
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    t.exports = r = i.extend({
        defaults: {
            searchBox: -10,
            overviewBox: 30,
            headerBox: -1,
            alignmentBody: 0,
            scaleSlider: 50
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    t.exports = r = i.extend({
        defaults: {
            sequences: !0,
            markers: !0,
            metacell: !1,
            conserv: !1,
            overviewbox: !1,
            seqlogo: !1,
            gapHeader: !1,
            leftHeader: !0,
            scaleslider: !1,
            labels: !0,
            labelName: !0,
            labelId: !0,
            labelPartition: !1,
            labelCheckbox: !1,
            metaGaps: !0,
            metaIdentity: !0,
            metaLinks: !0
        },
        constructor: function(t, e) {
            return this.calcDefaults(e.model), i.apply(this, arguments)
        },
        initialize: function() {
            return this.listenTo(this, "change:metaLinks change:metaIdentity change:metaGaps", function() {
                return this.trigger("change:metacell")
            }, this), this.listenTo(this, "change:labelName change:labelId change:labelPartition change:labelCheckbox", function() {
                return this.trigger("change:labels")
            }, this), this.listenTo(this, "change:markers change:conserv change:seqlogo change:gapHeader", function() {
                return this.trigger("change:header")
            }, this)
        },
        calcDefaults: function(t) {
            if (t.length > 0) {
                var e = t.at(0),
                    n = e.get("ids");
                if (void 0 !== n && 0 === Object.keys(n).length) return this.defaults.metaLinks = !1
            }
        }
    })
}, function(t, e, n) {
    "use strict";
    var r, i = n(1).Model;
    t.exports = r = i.extend({
        constructor: function(t, e) {
            return this.calcDefaults(e.model), i.apply(this, arguments), this.g = e.g, this.listenTo(this, "change:labelIdLength change:labelNameLength change:labelPartLength change:labelCheckLength", function() {
                return this.trigger("change:labelWidth", this.getLabelWidth())
            }, this), this.listenTo(this, "change:metaLinksWidth change:metaIdentWidth change:metaGapWidth", function() {
                return this.trigger("change:metaWidth", this.getMetaWidth())
            }, this), this
        },
        defaults: {
            alignmentWidth: "auto",
            alignmentHeight: 225,
            columnWidth: 15,
            rowHeight: 15,
            autoResize: !0,
            labelIdLength: 20,
            labelNameLength: 100,
            labelPartLength: 15,
            labelCheckLength: 15,
            labelFontsize: 13,
            labelLineHeight: "13px",
            markerFontsize: "10px",
            stepSize: 1,
            markerStepSize: 2,
            markerHeight: 20,
            residueFont: "13",
            canvasEventScale: 1,
            minLetterDrawSize: 11,
            boxRectHeight: 2,
            boxRectWidth: 2,
            overviewboxPaddingTop: 10,
            metaGapWidth: 35,
            metaIdentWidth: 40,
            _alignmentScrollLeft: 0,
            _alignmentScrollTop: 0
        },
        calcDefaults: function(t) {
            return t.getMaxLength() < 200 && t.length < 30 && (this.defaults.boxRectWidth = this.defaults.boxRectHeight = 5), this
        },
        getAlignmentWidth: function(t) {
            return this.get("autoResize") && void 0 !== t ? this.get("columnWidth") * t : void 0 === this.get("alignmentWidth") || "auto" === this.get("alignmentWidth") || 0 === this.get("alignmentWidth") ? this._adjustWidth() : this.get("alignmentWidth")
        },
        setLeftOffset: function(t) {
            var e = t;
            return e = Math.max(0, e), e -= this.g.columns.calcHiddenColumns(e), this.set("_alignmentScrollLeft", e * this.get("columnWidth"))
        },
        setTopOffset: function(t) {
            for (var e = Math.max(0, t - 1), n = 0, r = 0; 0 < e ? r <= e : r >= e; 0 < e ? r++ : r--) n += this.model.at(r).attributes.height || 1;
            return this.set("_alignmentScrollTop", n * this.get("rowHeight"))
        },
        getLeftBlockWidth: function() {
            var t = 0;
            return this.g.vis.get("labels") && (t += this.getLabelWidth()), this.g.vis.get("metacell") && (t += this.getMetaWidth()), t
        },
        getMetaWidth: function() {
            var t = 0;
            return this.g.vis.get("metaGaps") && (t += this.get("metaGapWidth")), this.g.vis.get("metaIdentity") && (t += this.get("metaIdentWidth")), this.g.vis.get("metaLinks") && (t += this.get("metaLinksWidth")), t
        },
        getLabelWidth: function() {
            var t = 0;
            return this.g.vis.get("labelName") && (t += this.get("labelNameLength")), this.g.vis.get("labelId") && (t += this.get("labelIdLength")), this.g.vis.get("labelPartition") && (t += this.get("labelPartLength")), this.g.vis.get("labelCheckbox") && (t += this.get("labelCheckLength")), t
        },
        _adjustWidth: function() {
            if (void 0 !== this.el && void 0 !== this.model) {
                if (null != this.el.parentNode && 0 !== this.el.parentNode.offsetWidth) var t = this.el.parentNode.offsetWidth;
                else t = document.body.clientWidth - 35;
                var e = t - this.getLeftBlockWidth(),
                    n = this.getAlignmentWidth(this.model.getMaxLength() - this.g.columns.get("hidden").length),
                    r = Math.min(e, n);
                return r = Math.floor(r / this.get("columnWidth")) * this.get("columnWidth"), this.attributes.alignmentWidth = r
            }
        },
        autoResize: function() {
            if (this.get("autoResize")) return this._adjustWidth(this.el, this.model)
        },
        autoHeight: function(t) {
            var e = this.getMaxAlignmentHeight();
            return void 0 !== t && t > 0 && (e = Math.min(e, t)), this.set("alignmentHeight", e)
        },
        setEl: function(t, e) {
            return this.el = t, this.model = e
        },
        _checkScrolling: function(t, e) {
            var n = t[0],
                r = t[1];
            return this.set("_alignmentScrollLeft", n, e), this.set("_alignmentScrollTop", r, e)
        },
        getMaxAlignmentHeight: function() {
            var t = 0;
            return this.model.each(function(e) {
                return t += e.attributes.height || 1
            }), t * this.get("rowHeight")
        },
        getMaxAlignmentWidth: function() {
            return this.model.getMaxLength() * this.get("columnWidth")
        }
    })
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(3),
        o = n(50),
        s = r(o),
        u = n(1).Collection,
        a = u.extend({
            model: s["default"],
            constructor: function() {
                return this.startOnCache = [], this.on("all", function() {
                    return this.startOnCache = []
                }, this), u.apply(this, arguments)
            },
            startOn: function(t) {
                return null == this.startOnCache[t] && (this.startOnCache[t] = this.where({
                    xStart: t
                })), this.startOnCache[t]
            },
            contains: function(t) {
                return this.reduce(function(e, n) {
                    return n || e.contains(t)
                }, !1)
            },
            getFeatureOnRow: function(t, e) {
                return this.filter(function(n) {
                    return n.get("row") === t && n.get("xStart") <= e && e <= n.get("xEnd")
                })
            },
            assignRows: function() {
                var t = this.max(function(t) {
                        return t.get("xEnd")
                    }).attributes.xEnd,
                    e = function() {
                        for (var e = [], n = 0; 0 < t ? n <= t : n >= t; 0 < t ? n++ : n--) e.push(0);
                        return e
                    }();
                return this.each(function(t) {
                    for (var n = 0, r = t.get("xStart"), i = t.get("xEnd"), o = r; r < i ? o <= i : o >= i; r < i ? o++ : o--) e[o] > n && (n = e[o]), e[o]++;
                    return t.set("row", n)
                }), (0, i.max)(e)
            },
            getCurrentHeight: function() {
                return this.max(function(t) {
                    return t.get("row")
                }).attributes.row + 1
            },
            getMinRows: function() {
                var t = this.max(function(t) {
                        return t.get("xEnd")
                    }).attributes.xEnd,
                    e = function() {
                        for (var e = [], n = 0; 0 < t ? n <= t : n >= t; 0 < t ? n++ : n--) e.push(0);
                        return e
                    }();
                return this.each(function(t) {
                    return function() {
                        for (var n = [], r = t.get("xStart"), i = t.get("xEnd"), o = r; r < i ? o <= i : o >= i; r < i ? o++ : o++) n.push(e[o]++);
                        return n
                    }()
                }), (0, i.max)(e)
            }
        });
    e["default"] = a
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(36),
        o = r(i),
        s = n(34),
        u = r(s),
        a = n(1).Collection,
        l = a.extend({
            model: o["default"],
            constructor: function(t, e) {
                var n = this;
                return a.apply(this, arguments), this.g = e, this.on("add reset remove", function() {
                    return n.lengthCache = null, n._bindSeqsWithFeatures()
                }, this), this.on("reset", function() {
                    return n._autoSetRefSeq()
                }), this._autoSetRefSeq(), this.lengthCache = null, this.features = {}, this
            },
            getMaxLength: function() {
                return 0 === this.models.length ? 0 : (null === this.lengthCache && (this.lengthCache = this.max(function(t) {
                    return t.get("seq").length
                }).get("seq").length), this.lengthCache)
            },
            prev: function(t, e) {
                var n = this.indexOf(t) - 1;
                return n < 0 && e && (n = this.length - 1), this.at(n)
            },
            next: function(t, e) {
                var n = this.indexOf(t) + 1;
                return n === this.length && e && (n = 0), this.at(n)
            },
            calcHiddenSeqs: function(t) {
                for (var e = t, n = 0; 0 < e ? n <= e : n >= e; 0 < e ? n++ : n--) this.at(n).get("hidden") && e++;
                return e - t
            },
            addFeatures: function(t) {
                var e = this;
                if (null != t.config) {
                    var n = t;
                    t = t.seqs, null != n.config.colors && ! function() {
                        var e = n.config.colors;
                        _.each(t, function(t) {
                            return _.each(t, function(t) {
                                if (null != e[t.feature]) return t.fillColor = e[t.feature]
                            })
                        })
                    }()
                }
                return _.isEmpty(this.features) ? this.features = t : _.each(t, function(t, n) {
                    return e.features.hasOwnProperty(n) ? e.features[n] = _.union(e.features[n], t) : e.features[n] = t
                }), this._bindSeqsWithFeatures()
            },
            _bindSeqWithFeatures: function(t) {
                var e = this.features[t.attributes.name];
                e && (t.attributes.features = new u["default"](e), t.attributes.features.assignRows(), t.attributes.height = t.attributes.features.getCurrentHeight() + 1)
            },
            _bindSeqsWithFeatures: function() {
                var t = this;
                return this.each(function(e) {
                    return t._bindSeqWithFeatures(e)
                })
            },
            removeAllFeatures: function() {
                return delete this.features
            },
            _autoSetRefSeq: function() {
                if (this.length > 0) return this.at(0).set("ref", !0)
            },
            setRef: function(t) {
                var e = this.get(t);
                return this.each(function(n) {
                    if (t.cid) return e.cid === n.cid ? n.set("ref", !0) : n.set("ref", !1)
                }), this.g.config.set("hasRef", !0), this.trigger("change:reference", t)
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(34),
        o = r(i),
        s = n(1).Model,
        u = s.extend({
            defaults: {
                name: "",
                id: "",
                seq: "",
                height: 1,
                ref: !1
            },
            initialize: function() {
                if (this.set("grey", []), null == this.get("features")) return this.set("features", new o["default"])
            }
        });
    e["default"] = u
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e.exporter = e.file = e.seqgen = e.proxy = e.bmath = void 0;
    var i = n(15),
        o = r(i),
        s = n(19),
        u = r(s),
        a = n(39),
        l = r(a),
        c = n(17),
        f = r(c),
        h = n(16),
        d = r(h);
    e.bmath = o["default"], e.proxy = u["default"], e.seqgen = l["default"], e.file = f["default"], e.exporter = d["default"]
}, function(t, e) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e["default"] = function(t, e) {
        for (var n = e.name || e.url || "", r = n.split("."), i = r[r.length - 1] || "", o = 0; o < s.length; o++) {
            var u = s[o](t, i);
            if (u) return u
        }
        return "unknown"
    };
    var n = function(t, e) {
            return ("CLUSTAL" === t.substring(0, 7) || "clustal" == e || "aln" == e) && "clustal"
        },
        r = function(t, e) {
            return (">" === t.substring(0, 1) || "fasta" == e || "fa" == e) && "fasta"
        },
        i = function(t, e) {
            return ("(" === t.substring(0, 1) || "nwk" == e) && "newick"
        },
        o = function(t, e) {
            if (t.length <= 10) return !1;
            var n = t.split("\n");
            return n[0].indexOf("gff") >= 0 || e.indexOf("gff") >= 0 ? "gff" : n[0].indexOf("#") < 0 && 2 === n[0].split("\t").length && "gff"
        },
        s = [n, r, i, o]
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(15),
        o = r(i),
        s = n(72).seq,
        u = n(48),
        a = {
            _generateSequence: function(t) {
                for (var e = "", n = t - 1, r = 0; 0 < n ? r <= n : r >= n; 0 < n ? r++ : r--) e += a.getRandomChar();
                return e
            },
            getDummySequences: function(t, e) {
                var n = [];
                "undefined" != typeof t && null !== t || (t = o["default"].getRandomInt(3, 5)), "undefined" != typeof e && null !== e || (e = o["default"].getRandomInt(50, 200));
                for (var r = 1; 1 < t ? r <= t : r >= t; 1 < t ? r++ : r--) n.push(new s(a._generateSequence(e), "seq" + r, "r" + r));
                return n
            },
            getRandomChar: function(t) {
                var e = t || "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
                return e.charAt(Math.floor(Math.random() * e.length))
            },
            genConservedSequences: function(t, e, n) {
                var r = [];
                "undefined" != typeof t && null !== t || (t = o["default"].getRandomInt(3, 5)), "undefined" != typeof e && null !== e || (e = o["default"].getRandomInt(50, 200)), n = n || "ACDEFGHIKLMNPQRSTVWY---";
                for (var i = 1; 1 < t ? i <= t : i >= t; 1 < t ? i++ : i--) r[i - 1] = "";
                for (var l = .2, c = 1, f = e - 1, h = 0; 0 < f ? h <= f : h >= f; 0 < f ? h++ : h--) {
                    h % 3 === 0 && (c = o["default"].getRandomInt(50, 100) / 100);
                    for (var d = [], p = t - 1, g = 0; 0 < p ? g <= p : g >= p; 0 < p ? g++ : g--) {
                        for (var v = 0, m = void 0; v < 100;) {
                            m = a.getRandomChar(n);
                            var y = u(d);
                            if (y.addSeq(m), v++, Math.abs(c - y.scale(y.conservation())[0]) < l) break
                        }
                        r[g] += m, d.push(m)
                    }
                }
                for (var _ = [], b = 1; 1 < t ? b <= t : b >= t; 1 < t ? b++ : b--) _.push(new s(r[b - 1], "seq" + b, "r" + b));
                return _
            }
        };
    e["default"] = a
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(3),
        o = n(35),
        s = (r(o), function(t) {
            return this.msa = t, this
        }),
        u = {
            loadTree: function(t) {
                return this.msa.g["package"].loadPackages(["msa-tnt", "biojs-io-newick"], t)
            },
            showTree: function(t) {
                function e(t) {
                    if (null != t.children) t.children.forEach(function(t) {
                        return e(t)
                    });
                    else {
                        var n = a.filter(function(e) {
                            return e.name === t.name
                        })[0];
                        null != n && ("number" == typeof n.id ? (n.ids = ["s" + (n.id + 1)], t.name = "s" + (n.id + 1)) : t.name = n.id)
                    }
                }
                var r = window.require("biojs-io-newick"),
                    i = window.require("msa-tnt");
                if ("string" == typeof t) var o = r.parse_newick(t);
                else o = t;
                var s, u = new i.selections;
                0 === this.msa.el.getElementsByClassName("tnt_groupDiv").length ? (s = document.createElement("div"), this.msa.el.appendChild(s)) : (s = this.msa.el.getElementsByClassName("tnt_groupDiv")[0].parentNode, s.innerHTML = "");
                var a = this.msa.seqs.toJSON();
                e(o);
                var l = i.app({
                    seqs: a,
                    tree: o
                });
                return new i.adapters.tree({
                    model: l,
                    el: s,
                    sel: u
                }), new i.adapters.msa({
                    model: l,
                    sel: u,
                    msa: this.msa
                }), l.models.forEach(function(t) {
                    return delete t.collection, Object.setPrototypeOf(t, n(1).Model.prototype)
                }), this.msa.seqs.reset(l.models), console.log(this.msa.seqs)
            },
            require: function(t) {
                function e(e) {
                    return t.apply(this, arguments)
                }
                return e.toString = function() {
                    return t.toString()
                }, e
            }(function(t) {
                return n(156)(t)
            })
        };
    (0, i.extend)(s.prototype, u), e["default"] = s
}, function(t, e, n) {
    "use strict";
    var r = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        },
        i = n(22),
        o = n(23),
        s = n(42),
        u = function(t, e) {
            var n = t || {};
            e || (e = {}), this.cid = s.uniqueId("c"), this.attributes = {}, e.collection && (this.collection = e.collection), e.parse && (n = this.parse(n, e) || {}), n = s.defaults({}, n, s.result(this, "defaults")), this.set(n, e), this.changed = {}, this.initialize.apply(this, arguments)
        };
    s.extend(u.prototype, i, {
        changed: null,
        validationError: null,
        idAttribute: "id",
        initialize: function() {},
        toJSON: function(t) {
            return s.clone(this.attributes)
        },
        sync: function() {
            return Backbone.sync.apply(this, arguments)
        },
        get: function(t) {
            return this.attributes[t]
        },
        escape: function(t) {
            return s.escape(this.get(t))
        },
        has: function(t) {
            return null != this.get(t)
        },
        set: function(t, e, n) {
            var i, o, u, a, l, c, f, h;
            if (null == t) return this;
            if ("object" === ("undefined" == typeof t ? "undefined" : r(t)) ? (o = t, n = e) : (o = {})[t] = e, n || (n = {}), !this._validate(o, n)) return !1;
            u = n.unset, l = n.silent, a = [], c = this._changing, this._changing = !0, c || (this._previousAttributes = s.clone(this.attributes), this.changed = {}), h = this.attributes, f = this._previousAttributes, this.idAttribute in o && (this.id = o[this.idAttribute]);
            for (i in o) e = o[i], s.isEqual(h[i], e) || a.push(i), s.isEqual(f[i], e) ? delete this.changed[i] : this.changed[i] = e, u ? delete h[i] : h[i] = e;
            if (!l) {
                a.length && (this._pending = n);
                for (var d = 0, p = a.length; d < p; d++) this.trigger("change:" + a[d], this, h[a[d]], n)
            }
            if (c) return this;
            if (!l)
                for (; this._pending;) n = this._pending, this._pending = !1, this.trigger("change", this, n);
            return this._pending = !1, this._changing = !1, this
        },
        unset: function(t, e) {
            return this.set(t, void 0, s.extend({}, e, {
                unset: !0
            }))
        },
        clear: function(t) {
            var e = {};
            for (var n in this.attributes) e[n] = void 0;
            return this.set(e, s.extend({}, t, {
                unset: !0
            }))
        },
        hasChanged: function(t) {
            return null == t ? !s.isEmpty(this.changed) : s.has(this.changed, t)
        },
        changedAttributes: function(t) {
            if (!t) return !!this.hasChanged() && s.clone(this.changed);
            var e, n = !1,
                r = this._changing ? this._previousAttributes : this.attributes;
            for (var i in t) s.isEqual(r[i], e = t[i]) || ((n || (n = {}))[i] = e);
            return n
        },
        previous: function(t) {
            return null != t && this._previousAttributes ? this._previousAttributes[t] : null
        },
        previousAttributes: function() {
            return s.clone(this._previousAttributes)
        },
        fetch: function(t) {
            t = t ? s.clone(t) : {}, void 0 === t.parse && (t.parse = !0);
            var e = this,
                n = t.success;
            return t.success = function(r) {
                return !!e.set(e.parse(r, t), t) && (n && n(e, r, t), void e.trigger("sync", e, r, t))
            }, wrapError(this, t), this.sync("read", this, t)
        },
        save: function(t, e, n) {
            var i, o, u, a = this.attributes;
            if (null == t || "object" === ("undefined" == typeof t ? "undefined" : r(t)) ? (i = t, n = e) : (i = {})[t] = e, n = s.extend({
                    validate: !0
                }, n), i && !n.wait) {
                if (!this.set(i, n)) return !1
            } else if (!this._validate(i, n)) return !1;
            i && n.wait && (this.attributes = s.extend({}, a, i)), void 0 === n.parse && (n.parse = !0);
            var l = this,
                c = n.success;
            return n.success = function(t) {
                l.attributes = a;
                var e = l.parse(t, n);
                return n.wait && (e = s.extend(i || {}, e)), !(s.isObject(e) && !l.set(e, n)) && (c && c(l, t, n), void l.trigger("sync", l, t, n))
            }, wrapError(this, n), o = this.isNew() ? "create" : n.patch ? "patch" : "update", "patch" !== o || n.attrs || (n.attrs = i), u = this.sync(o, this, n), i && n.wait && (this.attributes = a), u
        },
        destroy: function l(t) {
            t = t ? s.clone(t) : {};
            var e = this,
                n = t.success,
                l = function() {
                    e.stopListening(), e.trigger("destroy", e, e.collection, t)
                };
            if (t.success = function(r) {
                    (t.wait || e.isNew()) && l(), n && n(e, r, t), e.isNew() || e.trigger("sync", e, r, t)
                }, this.isNew()) return t.success(), !1;
            wrapError(this, t);
            var r = this.sync("delete", this, t);
            return t.wait || l(), r
        },
        url: function() {
            var t = s.result(this, "urlRoot") || s.result(this.collection, "url") || urlError();
            return this.isNew() ? t : t.replace(/([^\/])$/, "$1/") + encodeURIComponent(this.id)
        },
        parse: function(t, e) {
            return t
        },
        clone: function() {
            return new this.constructor(this.attributes)
        },
        isNew: function() {
            return !this.has(this.idAttribute)
        },
        isValid: function(t) {
            return this._validate({}, s.extend(t || {}, {
                validate: !0
            }))
        },
        _validate: function(t, e) {
            if (!e.validate || !this.validate) return !0;
            t = s.extend({}, this.attributes, t);
            var n = this.validationError = this.validate(t, e) || null;
            return !n || (this.trigger("invalid", this, n, s.extend(e, {
                validationError: n
            })), !1)
        }
    });
    var a = ["keys", "values", "pairs", "invert", "pick", "omit", "chain", "isEmpty"];
    s.each(a, function(t) {
        s[t] && (u.prototype[t] = function() {
            var e = slice.call(arguments);
            return e.unshift(this.attributes), s[t].apply(s, e)
        })
    }), u.extend = o, t.exports = u
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    (function() {
        function n(t) {
            function e(e, n, r, i, o, s) {
                for (; o >= 0 && o < s; o += t) {
                    var u = i ? i[o] : o;
                    r = n(r, e[u], u, e)
                }
                return r
            }
            return function(n, r, i, o) {
                r = k(r, o, 4);
                var s = !C(n) && S.keys(n),
                    u = (s || n).length,
                    a = t > 0 ? 0 : u - 1;
                return arguments.length < 3 && (i = n[s ? s[a] : a], a += t), e(n, r, i, s, a, u)
            }
        }

        function s(t) {
            return function(e, n, r) {
                n = j(n, r);
                for (var i = A(e), o = t > 0 ? 0 : i - 1; o >= 0 && o < i; o += t)
                    if (n(e[o], o, e)) return o;
                return -1
            }
        }

        function u(t, e, n) {
            return function(r, i, o) {
                var s = 0,
                    u = A(r);
                if ("number" == typeof o) t > 0 ? s = o >= 0 ? o : Math.max(o + u, s) : u = o >= 0 ? Math.min(o + 1, u) : o + u + 1;
                else if (n && o && u) return o = n(r, i), r[o] === i ? o : -1;
                if (i !== i) return o = e(g.call(r, s, u), S.isNaN), o >= 0 ? o + s : -1;
                for (o = t > 0 ? s : u - 1; o >= 0 && o < u; o += t)
                    if (r[o] === i) return o;
                return -1
            }
        }

        function a(t, e) {
            var n = R.length,
                r = t.constructor,
                i = S.isFunction(r) && r.prototype || h,
                o = "constructor";
            for (S.has(t, o) && !S.contains(e, o) && e.push(o); n--;) o = R[n], o in t && t[o] !== i[o] && !S.contains(e, o) && e.push(o)
        }
        var l = this,
            c = l._,
            f = Array.prototype,
            h = Object.prototype,
            d = Function.prototype,
            p = f.push,
            g = f.slice,
            v = h.toString,
            m = h.hasOwnProperty,
            y = Array.isArray,
            _ = Object.keys,
            b = d.bind,
            x = Object.create,
            w = function() {},
            S = function G(t) {
                return t instanceof G ? t : this instanceof G ? void(this._wrapped = t) : new G(t)
            };
        "undefined" != typeof t && t.exports && (e = t.exports = S), e._ = S, S.VERSION = "1.8.3";
        var k = function(t, e, n) {
                if (void 0 === e) return t;
                switch (null == n ? 3 : n) {
                    case 1:
                        return function(n) {
                            return t.call(e, n)
                        };
                    case 2:
                        return function(n, r) {
                            return t.call(e, n, r)
                        };
                    case 3:
                        return function(n, r, i) {
                            return t.call(e, n, r, i)
                        };
                    case 4:
                        return function(n, r, i, o) {
                            return t.call(e, n, r, i, o)
                        }
                }
                return function() {
                    return t.apply(e, arguments)
                }
            },
            j = function(t, e, n) {
                return null == t ? S.identity : S.isFunction(t) ? k(t, e, n) : S.isObject(t) ? S.matcher(t) : S.property(t)
            };
        S.iteratee = function(t, e) {
            return j(t, e, 1 / 0)
        };
        var O = function(t, e) {
                return function(n) {
                    var r = arguments.length;
                    if (r < 2 || null == n) return n;
                    for (var i = 1; i < r; i++)
                        for (var o = arguments[i], s = t(o), u = s.length, a = 0; a < u; a++) {
                            var l = s[a];
                            e && void 0 !== n[l] || (n[l] = o[l])
                        }
                    return n
                }
            },
            E = function(t) {
                if (!S.isObject(t)) return {};
                if (x) return x(t);
                w.prototype = t;
                var e = new w;
                return w.prototype = null, e
            },
            M = function(t) {
                return function(e) {
                    return null == e ? void 0 : e[t]
                }
            },
            z = Math.pow(2, 53) - 1,
            A = M("length"),
            C = function(t) {
                var e = A(t);
                return "number" == typeof e && e >= 0 && e <= z
            };
        S.each = S.forEach = function(t, e, n) {
            e = k(e, n);
            var r, i;
            if (C(t))
                for (r = 0, i = t.length; r < i; r++) e(t[r], r, t);
            else {
                var o = S.keys(t);
                for (r = 0, i = o.length; r < i; r++) e(t[o[r]], o[r], t)
            }
            return t
        }, S.map = S.collect = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = Array(i), s = 0; s < i; s++) {
                var u = r ? r[s] : s;
                o[s] = e(t[u], u, t)
            }
            return o
        }, S.reduce = S.foldl = S.inject = n(1), S.reduceRight = S.foldr = n(-1), S.find = S.detect = function(t, e, n) {
            var r;
            if (r = C(t) ? S.findIndex(t, e, n) : S.findKey(t, e, n), void 0 !== r && r !== -1) return t[r]
        }, S.filter = S.select = function(t, e, n) {
            var r = [];
            return e = j(e, n), S.each(t, function(t, n, i) {
                e(t, n, i) && r.push(t)
            }), r
        }, S.reject = function(t, e, n) {
            return S.filter(t, S.negate(j(e)), n)
        }, S.every = S.all = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (!e(t[s], s, t)) return !1
            }
            return !0
        }, S.some = S.any = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (e(t[s], s, t)) return !0
            }
            return !1
        }, S.contains = S.includes = S.include = function(t, e, n, r) {
            return C(t) || (t = S.values(t)), ("number" != typeof n || r) && (n = 0), S.indexOf(t, e, n) >= 0
        }, S.invoke = function(t, e) {
            var n = g.call(arguments, 2),
                r = S.isFunction(e);
            return S.map(t, function(t) {
                var i = r ? e : t[e];
                return null == i ? i : i.apply(t, n)
            })
        }, S.pluck = function(t, e) {
            return S.map(t, S.property(e))
        }, S.where = function(t, e) {
            return S.filter(t, S.matcher(e))
        }, S.findWhere = function(t, e) {
            return S.find(t, S.matcher(e))
        }, S.max = function(t, e, n) {
            var r, i, o = -(1 / 0),
                s = -(1 / 0);
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r > o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i > s || i === -(1 / 0) && o === -(1 / 0)) && (o = t, s = i)
            });
            return o
        }, S.min = function(t, e, n) {
            var r, i, o = 1 / 0,
                s = 1 / 0;
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r < o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i < s || i === 1 / 0 && o === 1 / 0) && (o = t, s = i)
            });
            return o
        }, S.shuffle = function(t) {
            for (var e, n = C(t) ? t : S.values(t), r = n.length, i = Array(r), o = 0; o < r; o++) e = S.random(0, o), e !== o && (i[o] = i[e]), i[e] = n[o];
            return i
        }, S.sample = function(t, e, n) {
            return null == e || n ? (C(t) || (t = S.values(t)), t[S.random(t.length - 1)]) : S.shuffle(t).slice(0, Math.max(0, e))
        }, S.sortBy = function(t, e, n) {
            return e = j(e, n), S.pluck(S.map(t, function(t, n, r) {
                return {
                    value: t,
                    index: n,
                    criteria: e(t, n, r)
                }
            }).sort(function(t, e) {
                var n = t.criteria,
                    r = e.criteria;
                if (n !== r) {
                    if (n > r || void 0 === n) return 1;
                    if (n < r || void 0 === r) return -1
                }
                return t.index - e.index
            }), "value")
        };
        var T = function(t) {
            return function(e, n, r) {
                var i = {};
                return n = j(n, r), S.each(e, function(r, o) {
                    t(i, r, n(r, o, e))
                }), i
            }
        };
        S.groupBy = T(function(t, e, n) {
            S.has(t, n) ? t[n].push(e) : t[n] = [e]
        }), S.indexBy = T(function(t, e, n) {
            t[n] = e
        }), S.countBy = T(function(t, e, n) {
            S.has(t, n) ? t[n]++ : t[n] = 1
        }), S.toArray = function(t) {
            return t ? S.isArray(t) ? g.call(t) : C(t) ? S.map(t, S.identity) : S.values(t) : []
        }, S.size = function(t) {
            return null == t ? 0 : C(t) ? t.length : S.keys(t).length
        }, S.partition = function(t, e, n) {
            e = j(e, n);
            var r = [],
                i = [];
            return S.each(t, function(t, n, o) {
                (e(t, n, o) ? r : i).push(t)
            }), [r, i]
        }, S.first = S.head = S.take = function(t, e, n) {
            if (null != t) return null == e || n ? t[0] : S.initial(t, t.length - e)
        }, S.initial = function(t, e, n) {
            return g.call(t, 0, Math.max(0, t.length - (null == e || n ? 1 : e)))
        }, S.last = function(t, e, n) {
            if (null != t) return null == e || n ? t[t.length - 1] : S.rest(t, Math.max(0, t.length - e))
        }, S.rest = S.tail = S.drop = function(t, e, n) {
            return g.call(t, null == e || n ? 1 : e)
        }, S.compact = function(t) {
            return S.filter(t, S.identity)
        };
        var I = function K(t, e, n, r) {
            for (var i = [], o = 0, s = r || 0, u = A(t); s < u; s++) {
                var a = t[s];
                if (C(a) && (S.isArray(a) || S.isArguments(a))) {
                    e || (a = K(a, e, n));
                    var l = 0,
                        c = a.length;
                    for (i.length += c; l < c;) i[o++] = a[l++]
                } else n || (i[o++] = a)
            }
            return i
        };
        S.flatten = function(t, e) {
            return I(t, e, !1)
        }, S.without = function(t) {
            return S.difference(t, g.call(arguments, 1))
        }, S.uniq = S.unique = function(t, e, n, r) {
            S.isBoolean(e) || (r = n, n = e, e = !1), null != n && (n = j(n, r));
            for (var i = [], o = [], s = 0, u = A(t); s < u; s++) {
                var a = t[s],
                    l = n ? n(a, s, t) : a;
                e ? (s && o === l || i.push(a), o = l) : n ? S.contains(o, l) || (o.push(l), i.push(a)) : S.contains(i, a) || i.push(a)
            }
            return i
        }, S.union = function() {
            return S.uniq(I(arguments, !0, !0))
        }, S.intersection = function(t) {
            for (var e = [], n = arguments.length, r = 0, i = A(t); r < i; r++) {
                var o = t[r];
                if (!S.contains(e, o)) {
                    for (var s = 1; s < n && S.contains(arguments[s], o); s++);
                    s === n && e.push(o)
                }
            }
            return e
        }, S.difference = function(t) {
            var e = I(arguments, !0, !0, 1);
            return S.filter(t, function(t) {
                return !S.contains(e, t)
            })
        }, S.zip = function() {
            return S.unzip(arguments)
        }, S.unzip = function(t) {
            for (var e = t && S.max(t, A).length || 0, n = Array(e), r = 0; r < e; r++) n[r] = S.pluck(t, r);
            return n
        }, S.object = function(t, e) {
            for (var n = {}, r = 0, i = A(t); r < i; r++) e ? n[t[r]] = e[r] : n[t[r][0]] = t[r][1];
            return n
        }, S.findIndex = s(1), S.findLastIndex = s(-1), S.sortedIndex = function(t, e, n, r) {
            n = j(n, r, 1);
            for (var i = n(e), o = 0, s = A(t); o < s;) {
                var u = Math.floor((o + s) / 2);
                n(t[u]) < i ? o = u + 1 : s = u
            }
            return o
        }, S.indexOf = u(1, S.findIndex, S.sortedIndex), S.lastIndexOf = u(-1, S.findLastIndex), S.range = function(t, e, n) {
            null == e && (e = t || 0, t = 0), n = n || 1;
            for (var r = Math.max(Math.ceil((e - t) / n), 0), i = Array(r), o = 0; o < r; o++, t += n) i[o] = t;
            return i
        };
        var N = function(t, e, n, r, i) {
            if (!(r instanceof e)) return t.apply(n, i);
            var o = E(t.prototype),
                s = t.apply(o, i);
            return S.isObject(s) ? s : o
        };
        S.bind = function(t, e) {
            if (b && t.bind === b) return b.apply(t, g.call(arguments, 1));
            if (!S.isFunction(t)) throw new TypeError("Bind must be called on a function");
            var n = g.call(arguments, 2);
            return function r() {
                return N(t, r, e, this, n.concat(g.call(arguments)))
            }
        }, S.partial = function(t) {
            var e = g.call(arguments, 1);
            return function n() {
                for (var r = 0, i = e.length, o = Array(i), s = 0; s < i; s++) o[s] = e[s] === S ? arguments[r++] : e[s];
                for (; r < arguments.length;) o.push(arguments[r++]);
                return N(t, n, this, this, o)
            }
        }, S.bindAll = function(t) {
            var e, n, r = arguments.length;
            if (r <= 1) throw new Error("bindAll must be passed function names");
            for (e = 1; e < r; e++) n = arguments[e], t[n] = S.bind(t[n], t);
            return t
        }, S.memoize = function(t, e) {
            var n = function r(n) {
                var i = r.cache,
                    o = "" + (e ? e.apply(this, arguments) : n);
                return S.has(i, o) || (i[o] = t.apply(this, arguments)), i[o]
            };
            return n.cache = {}, n
        }, S.delay = function(t, e) {
            var n = g.call(arguments, 2);
            return setTimeout(function() {
                return t.apply(null, n)
            }, e)
        }, S.defer = S.partial(S.delay, S, 1), S.throttle = function(t, e, n) {
            var r, i, o, s = null,
                u = 0;
            n || (n = {});
            var a = function() {
                u = n.leading === !1 ? 0 : S.now(), s = null, o = t.apply(r, i), s || (r = i = null)
            };
            return function() {
                var l = S.now();
                u || n.leading !== !1 || (u = l);
                var c = e - (l - u);
                return r = this, i = arguments, c <= 0 || c > e ? (s && (clearTimeout(s), s = null), u = l, o = t.apply(r, i), s || (r = i = null)) : s || n.trailing === !1 || (s = setTimeout(a, c)), o
            }
        }, S.debounce = function(t, e, n) {
            var r, i, o, s, u, a = function l() {
                var a = S.now() - s;
                a < e && a >= 0 ? r = setTimeout(l, e - a) : (r = null, n || (u = t.apply(o, i), r || (o = i = null)))
            };
            return function() {
                o = this, i = arguments, s = S.now();
                var l = n && !r;
                return r || (r = setTimeout(a, e)), l && (u = t.apply(o, i), o = i = null), u
            }
        }, S.wrap = function(t, e) {
            return S.partial(e, t)
        }, S.negate = function(t) {
            return function() {
                return !t.apply(this, arguments)
            }
        }, S.compose = function() {
            var t = arguments,
                e = t.length - 1;
            return function() {
                for (var n = e, r = t[e].apply(this, arguments); n--;) r = t[n].call(this, r);
                return r
            }
        }, S.after = function(t, e) {
            return function() {
                if (--t < 1) return e.apply(this, arguments)
            }
        }, S.before = function(t, e) {
            var n;
            return function() {
                return --t > 0 && (n = e.apply(this, arguments)), t <= 1 && (e = null), n
            }
        }, S.once = S.partial(S.before, 2);
        var L = !{
                toString: null
            }.propertyIsEnumerable("toString"),
            R = ["valueOf", "isPrototypeOf", "toString", "propertyIsEnumerable", "hasOwnProperty", "toLocaleString"];
        S.keys = function(t) {
            if (!S.isObject(t)) return [];
            if (_) return _(t);
            var e = [];
            for (var n in t) S.has(t, n) && e.push(n);
            return L && a(t, e), e
        }, S.allKeys = function(t) {
            if (!S.isObject(t)) return [];
            var e = [];
            for (var n in t) e.push(n);
            return L && a(t, e), e
        }, S.values = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = t[e[i]];
            return r
        }, S.mapObject = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = i.length, s = {}, u = 0; u < o; u++) r = i[u], s[r] = e(t[r], r, t);
            return s
        }, S.pairs = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = [e[i], t[e[i]]];
            return r
        }, S.invert = function(t) {
            for (var e = {}, n = S.keys(t), r = 0, i = n.length; r < i; r++) e[t[n[r]]] = n[r];
            return e
        }, S.functions = S.methods = function(t) {
            var e = [];
            for (var n in t) S.isFunction(t[n]) && e.push(n);
            return e.sort()
        }, S.extend = O(S.allKeys), S.extendOwn = S.assign = O(S.keys), S.findKey = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = 0, s = i.length; o < s; o++)
                if (r = i[o], e(t[r], r, t)) return r
        }, S.pick = function(t, e, n) {
            var r, i, o = {},
                s = t;
            if (null == s) return o;
            S.isFunction(e) ? (i = S.allKeys(s), r = k(e, n)) : (i = I(arguments, !1, !1, 1), r = function(t, e, n) {
                return e in n
            }, s = Object(s));
            for (var u = 0, a = i.length; u < a; u++) {
                var l = i[u],
                    c = s[l];
                r(c, l, s) && (o[l] = c)
            }
            return o
        }, S.omit = function(t, e, n) {
            if (S.isFunction(e)) e = S.negate(e);
            else {
                var r = S.map(I(arguments, !1, !1, 1), String);
                e = function(t, e) {
                    return !S.contains(r, e)
                }
            }
            return S.pick(t, e, n)
        }, S.defaults = O(S.allKeys, !0), S.create = function(t, e) {
            var n = E(t);
            return e && S.extendOwn(n, e), n
        }, S.clone = function(t) {
            return S.isObject(t) ? S.isArray(t) ? t.slice() : S.extend({}, t) : t
        }, S.tap = function(t, e) {
            return e(t), t
        }, S.isMatch = function(t, e) {
            var n = S.keys(e),
                r = n.length;
            if (null == t) return !r;
            for (var i = Object(t), o = 0; o < r; o++) {
                var s = n[o];
                if (e[s] !== i[s] || !(s in i)) return !1
            }
            return !0
        };
        var q = function X(t, e, n, r) {
            if (t === e) return 0 !== t || 1 / t === 1 / e;
            if (null == t || null == e) return t === e;
            t instanceof S && (t = t._wrapped), e instanceof S && (e = e._wrapped);
            var i = v.call(t);
            if (i !== v.call(e)) return !1;
            switch (i) {
                case "[object RegExp]":
                case "[object String]":
                    return "" + t == "" + e;
                case "[object Number]":
                    return +t !== +t ? +e !== +e : 0 === +t ? 1 / +t === 1 / e : +t === +e;
                case "[object Date]":
                case "[object Boolean]":
                    return +t === +e
            }
            var s = "[object Array]" === i;
            if (!s) {
                if ("object" != ("undefined" == typeof t ? "undefined" : o(t)) || "object" != ("undefined" == typeof e ? "undefined" : o(e))) return !1;
                var u = t.constructor,
                    a = e.constructor;
                if (u !== a && !(S.isFunction(u) && u instanceof u && S.isFunction(a) && a instanceof a) && "constructor" in t && "constructor" in e) return !1
            }
            n = n || [], r = r || [];
            for (var l = n.length; l--;)
                if (n[l] === t) return r[l] === e;
            if (n.push(t), r.push(e), s) {
                if (l = t.length, l !== e.length) return !1;
                for (; l--;)
                    if (!X(t[l], e[l], n, r)) return !1
            } else {
                var c, f = S.keys(t);
                if (l = f.length, S.keys(e).length !== l) return !1;
                for (; l--;)
                    if (c = f[l], !S.has(e, c) || !X(t[c], e[c], n, r)) return !1
            }
            return n.pop(), r.pop(), !0
        };
        S.isEqual = function(t, e) {
            return q(t, e)
        }, S.isEmpty = function(t) {
            return null == t || (C(t) && (S.isArray(t) || S.isString(t) || S.isArguments(t)) ? 0 === t.length : 0 === S.keys(t).length)
        }, S.isElement = function(t) {
            return !(!t || 1 !== t.nodeType)
        }, S.isArray = y || function(t) {
            return "[object Array]" === v.call(t)
        }, S.isObject = function(t) {
            var e = "undefined" == typeof t ? "undefined" : o(t);
            return "function" === e || "object" === e && !!t
        }, S.each(["Arguments", "Function", "String", "Number", "Date", "RegExp", "Error"], function(t) {
            S["is" + t] = function(e) {
                return v.call(e) === "[object " + t + "]"
            }
        }), S.isArguments(arguments) || (S.isArguments = function(t) {
            return S.has(t, "callee")
        }), "function" != typeof /./ && "object" != ("undefined" == typeof Int8Array ? "undefined" : o(Int8Array)) && (S.isFunction = function(t) {
            return "function" == typeof t || !1
        }), S.isFinite = function(t) {
            return isFinite(t) && !isNaN(parseFloat(t))
        }, S.isNaN = function(t) {
            return S.isNumber(t) && t !== +t
        }, S.isBoolean = function(t) {
            return t === !0 || t === !1 || "[object Boolean]" === v.call(t)
        }, S.isNull = function(t) {
            return null === t
        }, S.isUndefined = function(t) {
            return void 0 === t
        }, S.has = function(t, e) {
            return null != t && m.call(t, e)
        }, S.noConflict = function() {
            return l._ = c, this
        }, S.identity = function(t) {
            return t
        }, S.constant = function(t) {
            return function() {
                return t
            }
        }, S.noop = function() {}, S.property = M, S.propertyOf = function(t) {
            return null == t ? function() {} : function(e) {
                return t[e]
            }
        }, S.matcher = S.matches = function(t) {
            return t = S.extendOwn({}, t),
                function(e) {
                    return S.isMatch(e, t)
                }
        }, S.times = function(t, e, n) {
            var r = Array(Math.max(0, t));
            e = k(e, n, 1);
            for (var i = 0; i < t; i++) r[i] = e(i);
            return r
        }, S.random = function(t, e) {
            return null == e && (e = t, t = 0), t + Math.floor(Math.random() * (e - t + 1))
        }, S.now = Date.now || function() {
            return (new Date).getTime()
        };
        var F = {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': "&quot;",
                "'": "&#x27;",
                "`": "&#x60;"
            },
            P = S.invert(F),
            B = function(t) {
                var e = function(e) {
                        return t[e]
                    },
                    n = "(?:" + S.keys(t).join("|") + ")",
                    r = RegExp(n),
                    i = RegExp(n, "g");
                return function(t) {
                    return t = null == t ? "" : "" + t, r.test(t) ? t.replace(i, e) : t
                }
            };
        S.escape = B(F), S.unescape = B(P), S.result = function(t, e, n) {
            var r = null == t ? void 0 : t[e];
            return void 0 === r && (r = n), S.isFunction(r) ? r.call(t) : r
        };
        var W = 0;
        S.uniqueId = function(t) {
            var e = ++W + "";
            return t ? t + e : e
        }, S.templateSettings = {
            evaluate: /<%([\s\S]+?)%>/g,
            interpolate: /<%=([\s\S]+?)%>/g,
            escape: /<%-([\s\S]+?)%>/g
        };
        var D = /(.)^/,
            H = {
                "'": "'",
                "\\": "\\",
                "\r": "r",
                "\n": "n",
                "\u2028": "u2028",
                "\u2029": "u2029"
            },
            U = /\\|'|\r|\n|\u2028|\u2029/g,
            V = function(t) {
                return "\\" + H[t]
            };
        S.template = function(t, e, n) {
            !e && n && (e = n), e = S.defaults({}, e, S.templateSettings);
            var r = RegExp([(e.escape || D).source, (e.interpolate || D).source, (e.evaluate || D).source].join("|") + "|$", "g"),
                i = 0,
                o = "__p+='";
            t.replace(r, function(e, n, r, s, u) {
                return o += t.slice(i, u).replace(U, V), i = u + e.length, n ? o += "'+\n((__t=(" + n + "))==null?'':_.escape(__t))+\n'" : r ? o += "'+\n((__t=(" + r + "))==null?'':__t)+\n'" : s && (o += "';\n" + s + "\n__p+='"),
                    e
            }), o += "';\n", e.variable || (o = "with(obj||{}){\n" + o + "}\n"), o = "var __t,__p='',__j=Array.prototype.join,print=function(){__p+=__j.call(arguments,'');};\n" + o + "return __p;\n";
            try {
                var s = new Function(e.variable || "obj", "_", o)
            } catch (u) {
                throw u.source = o, u
            }
            var a = function(t) {
                    return s.call(this, t, S)
                },
                l = e.variable || "obj";
            return a.source = "function(" + l + "){\n" + o + "}", a
        }, S.chain = function(t) {
            var e = S(t);
            return e._chain = !0, e
        };
        var $ = function(t, e) {
            return t._chain ? S(e).chain() : e
        };
        S.mixin = function(t) {
            S.each(S.functions(t), function(e) {
                var n = S[e] = t[e];
                S.prototype[e] = function() {
                    var t = [this._wrapped];
                    return p.apply(t, arguments), $(this, n.apply(S, t))
                }
            })
        }, S.mixin(S), S.each(["pop", "push", "reverse", "shift", "sort", "splice", "unshift"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                var n = this._wrapped;
                return e.apply(n, arguments), "shift" !== t && "splice" !== t || 0 !== n.length || delete n[0], $(this, n)
            }
        }), S.each(["concat", "join", "slice"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                return $(this, e.apply(this._wrapped, arguments))
            }
        }), S.prototype.value = function() {
            return this._wrapped
        }, S.prototype.valueOf = S.prototype.toJSON = S.prototype.value, S.prototype.toString = function() {
            return "" + this._wrapped
        }, r = [], i = function() {
            return S
        }.apply(e, r), !(void 0 !== i && (t.exports = i))
    }).call(void 0)
}, function(t, e) {
    "use strict";

    function n(t) {
        var e = {};
        return t.split(";").forEach(function(t) {
            var n, r, i;
            t.indexOf("=") > 0 ? (n = t.split("="), r = n[0], i = n[1], e[r] = i) : t.indexOf(" ") > 0 && (n = t.split(" "), r = n[0], i = n[1].replace(/"/g, ""), e[r] = i)
        }), e
    }

    function r(t) {
        var e = t.toString(16);
        return 1 === e.length ? "0" + e : e
    }

    function i(t, e, n) {
        return 3 === t.length ? i(t[0], t[1], t[2]) : "#" + r(t) + r(e) + r(n)
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e.extractKeys = n, e.rgbToHex = i, e["default"] = {
        extractKeys: n,
        rgbToHex: i
    }
}, function(t, e) {
    "use strict";
    t.exports = function(t, e, n) {
        t.beginPath(), t.moveTo(0, e), t.lineTo(n, e), t.lineWidth = 1, t.strokeStyle = "#999999", t.stroke()
    }
}, function(t, e) {
    "use strict";
    t.exports = function(t, e, n, r, i, o, s) {
        t.font = o + "px Arial", t.textAlign = s ? "right" : "center", t.fillStyle = "#666666", t.fillText(i, e + r / 2, n)
    }
}, function(t, e) {
    "use strict";
    t.exports = function(t, e, n, r, i) {
        i = i || "#999999", t.beginPath(), t.moveTo(e, n), t.lineTo(e, n + r), t.lineWidth = 1, t.strokeStyle = i, t.stroke()
    }
}, function(t, e) {
    "use strict";
    var n;
    t.exports = n = {
        rel: function(t) {
            var e, n, r, i;
            return e = t.offsetX, n = t.offsetY, void 0 == e && (r = i.getBoundingClientRect(), i = t.target || t.srcElement, void 0 == e && (e = t.clientX - r.left, n = t.clientY - r.top), void 0 == e && (e = t.pageX - i.offsetLeft, n = t.pageY - i.offsetTop), void 0 == e) ? void console.log(t, "no mouse event defined. your browser sucks") : [e, n]
        },
        abs: function(t) {
            var e, n;
            return e = t.pageX, n = t.pageY, void 0 == e && (e = t.layerX, n = t.layerY), void 0 == e && (e = t.clientX, n = t.clientY), void 0 == e && (e = t.x, n = t.y), [e, n]
        },
        wheelDelta: function(t) {
            var e;
            return e = [t.deltaX, t.deltaY], void 0 == e[0] && t.mozMovementX && (e = [0, t.mozMovementX]), isNaN(e[0]) && (e[0] = 0), isNaN(e[1]) && (e[1] = 0), e
        }
    }
}, function(t, e, n) {
    "use strict";
    var r = n(110),
        i = function s(t, e) {
            if (!this || this.constructor !== s) return new s(t);
            if (void 0 === t || "string" == typeof t) throw new TypeError("you need to give the seq stat an array");
            this.resetSeqs(t), this.alphabetSize = 4, this._useBackground = !1, this.useGaps = !1, this.ignoredChars = ["-", "*"], r.extend(this, e)
        };
    i.prototype.addSeq = function(t) {
        this.seqs.push(t), this._reset()
    }, i.prototype.removeSeq = function(t) {
        "number" == typeof t ? this.seqs.splice(t, 1) : r.each(this.seqs, function(e, n) {
            t === e && this.seqs.splice(n, 1)
        }.bind(this)), this._reset()
    }, i.prototype.addSeqs = function(t) {
        t.forEach(function(t) {
            this.addSeq(t)
        }.bind(this))
    }, i.prototype.resetSeqs = function(t) {
        if (this.seqs = [], !t instanceof Array || "at" in t) {
            this.mseqs = t;
            var e = function() {
                var t = this.mseqs.pluck("seq");
                this.resetSeqs(t)
            };
            t.on("add change reset ", e, this), e.call(this)
        } else this.addSeqs(t), this._reset(), this.trigger("reset")
    };
    var o = ["consensus", "frequency", "maxLength", "ic", "gaps"];
    i.prototype._reset = function() {
        for (var t = 0; t < o.length; t++) this["_" + o[t]] = void 0;
        this._identity = void 0, this._background = void 0
    }, i.prototype.setBackground = function(t) {
        this._useBackground = t, this._reset()
    }, i.prototype.useBackground = function() {
        this.setBackground(!0)
    }, i.prototype.setDNA = function() {
        this.alphabetSize = 4
    }, i.prototype.setProtein = function() {
        this.alphabetSize = 20
    }, o.forEach(function(t) {
        i.prototype[t] = function() {
            return void 0 === this["_" + t] && (this["_" + t] = this[t + "Calc"]()), this["_" + t]
        }
    }), i.prototype.identity = function(t) {
        var e;
        return (void 0 === this._identity || t) && (e = this.identityCalc(t), this._identity = void 0), this._identity || e
    }, i.prototype.background = function() {
        return void 0 !== this.bg ? this.bg : (void 0 === this._background && this.backgroundCalc(), this._background)
    }, i.prototype.frequencyCalc = function(t) {
        var e, n;
        e = new Array(this.maxLength()), n = new Array(this.seqs.length);
        var i = this.ignoredChars;
        return void 0 !== t && t.all && (i = []), r.each(this.seqs, function(t) {
            r.each(t, function(t, r) {
                i.indexOf(t) >= 0 || (void 0 === e[r] && (e[r] = {}), void 0 === e[r][t] && (e[r][t] = 0), e[r][t]++, void 0 === n[r] && (n[r] = 0), n[r]++)
            })
        }), r.each(e, function(t, i) {
            return r.each(t, function(t, r) {
                return e[i][r] = t / n[i]
            })
        }), this._frequency = e, e
    }, i.prototype.backgroundCalc = function() {
        var t = {},
            e = 0;
        return r.each(this.seqs, function(n) {
            r.each(n, function(n) {
                return void 0 === t[n] && (t[n] = 0), t[n]++, e++
            })
        }), t = r.mapValues(t, function(t) {
            return t / e
        }), this._background = t, t
    }, i.prototype.icCalc = function() {
        var t = this.frequency();
        if (this._useBackground) var e = this.background();
        var n = this.ignoredChars,
            i = this._useBackground,
            o = r.map(t, function(t) {
                return r.reduce(t, function(t, r, o) {
                    return n.indexOf(o) >= 0 ? t : (i && (r /= e[o]), t - r * (Math.log(r) / Math.log(2)))
                }, 0)
            });
        return this._ic = o, o
    }, i.prototype.conservation = function(t) {
        var e = this.ic(),
            n = this.gaps(),
            i = this;
        t = t || this.alphabetSize;
        var o = Math.log(t) / Math.log(2),
            s = 0;
        return r.map(e, function(t) {
            var e = o - t;
            return i.useGaps && (e *= 1 - n[s++]), e
        })
    }, i.prototype.conservResidue = function(t) {
        var e, n = t ? t.alphabetSize : void 0,
            i = this.ignoredChars;
        e = void 0 !== t && t.scaled ? this.scale(this.conservation(n)) : this.conservation(n);
        var o, s = this.frequency();
        return r.map(s, function(t, n) {
            o = r.reject(r.keys(t), function(t) {
                return i.indexOf(t) >= 0
            });
            var s = {};
            return r.each(o, function(r) {
                s[r] = t[r] * e[n]
            }), s
        })
    }, i.prototype.conservResidue2 = function(t) {
        var e = this.frequency(),
            n = this.conservation(t),
            i = this.background();
        return r.map(e, function(t, o) {
            return r.map(t, function(t) {
                var s = r.reduce(e[o], function(t, e) {
                    return t + e / i[o]
                }, 0);
                return t / i[o] / s * n[o]
            }, 0)
        })
    }, i.prototype.scale = function(t, e) {
        e = e || this.alphabetSize;
        var n = Math.log(e) / Math.log(2);
        return r.map(t, function(t) {
            return t / n
        })
    }, i.prototype.maxLengthCalc = function() {
        return 0 === this.seqs.length ? 0 : r.max(this.seqs, function(t) {
            return t.length
        }).length
    }, i.prototype.consensusCalc = function() {
        var t = new Array(this.maxLength());
        return r.each(this.seqs, function(e) {
            r.each(e, function(e, n) {
                void 0 === t[n] && (t[n] = {}), void 0 === t[n][e] && (t[n][e] = 0), t[n][e]++
            })
        }), this._consensus = r.reduce(t, function(t, e) {
            var n;
            return n = r.keys(e), t += r.max(n, function(t) {
                return e[t]
            })
        }, ""), this._consensus
    }, i.prototype.identityCalc = function(t) {
        var e = t || this.consensus();
        return this._identity = this.seqs.map(function(t) {
            for (var n = 0, r = 0, i = 0; i < t.length; i++) "-" !== t[i] && "-" !== e[i] && (r++, t[i] === e[i] && n++);
            return n / r
        }), this._identity
    }, i.prototype.gapsCalc = function() {
        var t = this.maxLength();
        if (t <= 1 || "undefined" == typeof t) return [];
        var e = new Array(this.maxLength());
        return r.each(this.seqs, function(t) {
            r.each(t, function(t, n) {
                void 0 === e[n] && (e[n] = {
                    g: 0,
                    t: 0
                }), t = "-" === t ? "g" : "t", e[n][t]++
            })
        }), this._gaps = r.map(e, function(t) {
            return t.g / (t.g + t.t)
        }), this._gaps
    }, r.mixin({
        mapValues: function(t, e) {
            return r.object(r.keys(t), r.map(t, e))
        }
    }), n(13).mixin(i.prototype), t.exports = i
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(9);
    Object.defineProperty(e, "SelectionManager", {
        enumerable: !0,
        get: function() {
            return r.SelectionManager
        }
    })
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(1).Model,
        i = r.extend({
            defaults: {
                xStart: -1,
                xEnd: -1,
                height: -1,
                text: "",
                fillColor: "red",
                fillOpacity: .5,
                type: "rectangle",
                borderSize: 1,
                borderColor: "black",
                borderOpacity: .5,
                validate: !0,
                row: 0
            },
            initialize: function(t) {
                if (null != t.start && this.set("xStart", t.start - 1), null != t.end && this.set("xEnd", t.end - 1), null != t.attributes && (null != t.attributes.Name && this.set("text", t.attributes.Name), null != t.attributes.Color && this.set("fillColor", t.attributes.Color)), this.attributes.xEnd < this.attributes.xStart && console.warn("invalid feature range for", this.attributes), !_.isNumber(this.attributes.xStart) || !_.isNumber(this.attributes.xEnd)) return console.warn("please provide numeric feature ranges", t), this.set("xStart", parseInt(this.attributes.xStart)), this.set("xEnd", parseInt(this.attributes.xEnd))
            },
            validate: function() {
                if (isNaN(this.attributes.xStart || isNaN(this.attributes.xEnd))) return "features need integer start and end."
            },
            contains: function(t) {
                return this.attributes.xStart <= t && t <= this.attributes.xEnd
            }
        });
    e["default"] = i
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(36);
    Object.defineProperty(e, "seq", {
        enumerable: !0,
        get: function() {
            return r(i)["default"]
        }
    });
    var o = n(35);
    Object.defineProperty(e, "seqcol", {
        enumerable: !0,
        get: function() {
            return r(o)["default"]
        }
    });
    var s = n(50);
    Object.defineProperty(e, "feature", {
        enumerable: !0,
        get: function() {
            return r(s)["default"]
        }
    });
    var u = n(34);
    Object.defineProperty(e, "featurecol", {
        enumerable: !0,
        get: function() {
            return r(u)["default"]
        }
    })
}, function(t, e) {
    "use strict";

    function n(t) {
        var e = r.call(t);
        return "[object Function]" === e || "function" == typeof t && "[object RegExp]" !== e || "undefined" != typeof window && (t === window.setTimeout || t === window.alert || t === window.confirm || t === window.prompt)
    }
    t.exports = n;
    var r = Object.prototype.toString
}, function(t, e) {
    (function(e) {
        t.exports = e
    }).call(e, {})
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    (function() {
        function n(t) {
            function e(e, n, r, i, o, s) {
                for (; o >= 0 && o < s; o += t) {
                    var u = i ? i[o] : o;
                    r = n(r, e[u], u, e)
                }
                return r
            }
            return function(n, r, i, o) {
                r = k(r, o, 4);
                var s = !C(n) && S.keys(n),
                    u = (s || n).length,
                    a = t > 0 ? 0 : u - 1;
                return arguments.length < 3 && (i = n[s ? s[a] : a], a += t), e(n, r, i, s, a, u)
            }
        }

        function s(t) {
            return function(e, n, r) {
                n = j(n, r);
                for (var i = A(e), o = t > 0 ? 0 : i - 1; o >= 0 && o < i; o += t)
                    if (n(e[o], o, e)) return o;
                return -1
            }
        }

        function u(t, e, n) {
            return function(r, i, o) {
                var s = 0,
                    u = A(r);
                if ("number" == typeof o) t > 0 ? s = o >= 0 ? o : Math.max(o + u, s) : u = o >= 0 ? Math.min(o + 1, u) : o + u + 1;
                else if (n && o && u) return o = n(r, i), r[o] === i ? o : -1;
                if (i !== i) return o = e(g.call(r, s, u), S.isNaN), o >= 0 ? o + s : -1;
                for (o = t > 0 ? s : u - 1; o >= 0 && o < u; o += t)
                    if (r[o] === i) return o;
                return -1
            }
        }

        function a(t, e) {
            var n = R.length,
                r = t.constructor,
                i = S.isFunction(r) && r.prototype || h,
                o = "constructor";
            for (S.has(t, o) && !S.contains(e, o) && e.push(o); n--;) o = R[n], o in t && t[o] !== i[o] && !S.contains(e, o) && e.push(o)
        }
        var l = this,
            c = l._,
            f = Array.prototype,
            h = Object.prototype,
            d = Function.prototype,
            p = f.push,
            g = f.slice,
            v = h.toString,
            m = h.hasOwnProperty,
            y = Array.isArray,
            _ = Object.keys,
            b = d.bind,
            x = Object.create,
            w = function() {},
            S = function G(t) {
                return t instanceof G ? t : this instanceof G ? void(this._wrapped = t) : new G(t)
            };
        "undefined" != typeof t && t.exports && (e = t.exports = S), e._ = S, S.VERSION = "1.8.3";
        var k = function(t, e, n) {
                if (void 0 === e) return t;
                switch (null == n ? 3 : n) {
                    case 1:
                        return function(n) {
                            return t.call(e, n)
                        };
                    case 2:
                        return function(n, r) {
                            return t.call(e, n, r)
                        };
                    case 3:
                        return function(n, r, i) {
                            return t.call(e, n, r, i)
                        };
                    case 4:
                        return function(n, r, i, o) {
                            return t.call(e, n, r, i, o)
                        }
                }
                return function() {
                    return t.apply(e, arguments)
                }
            },
            j = function(t, e, n) {
                return null == t ? S.identity : S.isFunction(t) ? k(t, e, n) : S.isObject(t) ? S.matcher(t) : S.property(t)
            };
        S.iteratee = function(t, e) {
            return j(t, e, 1 / 0)
        };
        var O = function(t, e) {
                return function(n) {
                    var r = arguments.length;
                    if (r < 2 || null == n) return n;
                    for (var i = 1; i < r; i++)
                        for (var o = arguments[i], s = t(o), u = s.length, a = 0; a < u; a++) {
                            var l = s[a];
                            e && void 0 !== n[l] || (n[l] = o[l])
                        }
                    return n
                }
            },
            E = function(t) {
                if (!S.isObject(t)) return {};
                if (x) return x(t);
                w.prototype = t;
                var e = new w;
                return w.prototype = null, e
            },
            M = function(t) {
                return function(e) {
                    return null == e ? void 0 : e[t]
                }
            },
            z = Math.pow(2, 53) - 1,
            A = M("length"),
            C = function(t) {
                var e = A(t);
                return "number" == typeof e && e >= 0 && e <= z
            };
        S.each = S.forEach = function(t, e, n) {
            e = k(e, n);
            var r, i;
            if (C(t))
                for (r = 0, i = t.length; r < i; r++) e(t[r], r, t);
            else {
                var o = S.keys(t);
                for (r = 0, i = o.length; r < i; r++) e(t[o[r]], o[r], t)
            }
            return t
        }, S.map = S.collect = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = Array(i), s = 0; s < i; s++) {
                var u = r ? r[s] : s;
                o[s] = e(t[u], u, t)
            }
            return o
        }, S.reduce = S.foldl = S.inject = n(1), S.reduceRight = S.foldr = n(-1), S.find = S.detect = function(t, e, n) {
            var r;
            if (r = C(t) ? S.findIndex(t, e, n) : S.findKey(t, e, n), void 0 !== r && r !== -1) return t[r]
        }, S.filter = S.select = function(t, e, n) {
            var r = [];
            return e = j(e, n), S.each(t, function(t, n, i) {
                e(t, n, i) && r.push(t)
            }), r
        }, S.reject = function(t, e, n) {
            return S.filter(t, S.negate(j(e)), n)
        }, S.every = S.all = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (!e(t[s], s, t)) return !1
            }
            return !0
        }, S.some = S.any = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (e(t[s], s, t)) return !0
            }
            return !1
        }, S.contains = S.includes = S.include = function(t, e, n, r) {
            return C(t) || (t = S.values(t)), ("number" != typeof n || r) && (n = 0), S.indexOf(t, e, n) >= 0
        }, S.invoke = function(t, e) {
            var n = g.call(arguments, 2),
                r = S.isFunction(e);
            return S.map(t, function(t) {
                var i = r ? e : t[e];
                return null == i ? i : i.apply(t, n)
            })
        }, S.pluck = function(t, e) {
            return S.map(t, S.property(e))
        }, S.where = function(t, e) {
            return S.filter(t, S.matcher(e))
        }, S.findWhere = function(t, e) {
            return S.find(t, S.matcher(e))
        }, S.max = function(t, e, n) {
            var r, i, o = -(1 / 0),
                s = -(1 / 0);
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r > o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i > s || i === -(1 / 0) && o === -(1 / 0)) && (o = t, s = i)
            });
            return o
        }, S.min = function(t, e, n) {
            var r, i, o = 1 / 0,
                s = 1 / 0;
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r < o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i < s || i === 1 / 0 && o === 1 / 0) && (o = t, s = i)
            });
            return o
        }, S.shuffle = function(t) {
            for (var e, n = C(t) ? t : S.values(t), r = n.length, i = Array(r), o = 0; o < r; o++) e = S.random(0, o), e !== o && (i[o] = i[e]), i[e] = n[o];
            return i
        }, S.sample = function(t, e, n) {
            return null == e || n ? (C(t) || (t = S.values(t)), t[S.random(t.length - 1)]) : S.shuffle(t).slice(0, Math.max(0, e))
        }, S.sortBy = function(t, e, n) {
            return e = j(e, n), S.pluck(S.map(t, function(t, n, r) {
                return {
                    value: t,
                    index: n,
                    criteria: e(t, n, r)
                }
            }).sort(function(t, e) {
                var n = t.criteria,
                    r = e.criteria;
                if (n !== r) {
                    if (n > r || void 0 === n) return 1;
                    if (n < r || void 0 === r) return -1
                }
                return t.index - e.index
            }), "value")
        };
        var T = function(t) {
            return function(e, n, r) {
                var i = {};
                return n = j(n, r), S.each(e, function(r, o) {
                    t(i, r, n(r, o, e))
                }), i
            }
        };
        S.groupBy = T(function(t, e, n) {
            S.has(t, n) ? t[n].push(e) : t[n] = [e]
        }), S.indexBy = T(function(t, e, n) {
            t[n] = e
        }), S.countBy = T(function(t, e, n) {
            S.has(t, n) ? t[n]++ : t[n] = 1
        }), S.toArray = function(t) {
            return t ? S.isArray(t) ? g.call(t) : C(t) ? S.map(t, S.identity) : S.values(t) : []
        }, S.size = function(t) {
            return null == t ? 0 : C(t) ? t.length : S.keys(t).length
        }, S.partition = function(t, e, n) {
            e = j(e, n);
            var r = [],
                i = [];
            return S.each(t, function(t, n, o) {
                (e(t, n, o) ? r : i).push(t)
            }), [r, i]
        }, S.first = S.head = S.take = function(t, e, n) {
            if (null != t) return null == e || n ? t[0] : S.initial(t, t.length - e)
        }, S.initial = function(t, e, n) {
            return g.call(t, 0, Math.max(0, t.length - (null == e || n ? 1 : e)))
        }, S.last = function(t, e, n) {
            if (null != t) return null == e || n ? t[t.length - 1] : S.rest(t, Math.max(0, t.length - e))
        }, S.rest = S.tail = S.drop = function(t, e, n) {
            return g.call(t, null == e || n ? 1 : e)
        }, S.compact = function(t) {
            return S.filter(t, S.identity)
        };
        var I = function K(t, e, n, r) {
            for (var i = [], o = 0, s = r || 0, u = A(t); s < u; s++) {
                var a = t[s];
                if (C(a) && (S.isArray(a) || S.isArguments(a))) {
                    e || (a = K(a, e, n));
                    var l = 0,
                        c = a.length;
                    for (i.length += c; l < c;) i[o++] = a[l++]
                } else n || (i[o++] = a)
            }
            return i
        };
        S.flatten = function(t, e) {
            return I(t, e, !1)
        }, S.without = function(t) {
            return S.difference(t, g.call(arguments, 1))
        }, S.uniq = S.unique = function(t, e, n, r) {
            S.isBoolean(e) || (r = n, n = e, e = !1), null != n && (n = j(n, r));
            for (var i = [], o = [], s = 0, u = A(t); s < u; s++) {
                var a = t[s],
                    l = n ? n(a, s, t) : a;
                e ? (s && o === l || i.push(a), o = l) : n ? S.contains(o, l) || (o.push(l), i.push(a)) : S.contains(i, a) || i.push(a)
            }
            return i
        }, S.union = function() {
            return S.uniq(I(arguments, !0, !0))
        }, S.intersection = function(t) {
            for (var e = [], n = arguments.length, r = 0, i = A(t); r < i; r++) {
                var o = t[r];
                if (!S.contains(e, o)) {
                    for (var s = 1; s < n && S.contains(arguments[s], o); s++);
                    s === n && e.push(o)
                }
            }
            return e
        }, S.difference = function(t) {
            var e = I(arguments, !0, !0, 1);
            return S.filter(t, function(t) {
                return !S.contains(e, t)
            })
        }, S.zip = function() {
            return S.unzip(arguments)
        }, S.unzip = function(t) {
            for (var e = t && S.max(t, A).length || 0, n = Array(e), r = 0; r < e; r++) n[r] = S.pluck(t, r);
            return n
        }, S.object = function(t, e) {
            for (var n = {}, r = 0, i = A(t); r < i; r++) e ? n[t[r]] = e[r] : n[t[r][0]] = t[r][1];
            return n
        }, S.findIndex = s(1), S.findLastIndex = s(-1), S.sortedIndex = function(t, e, n, r) {
            n = j(n, r, 1);
            for (var i = n(e), o = 0, s = A(t); o < s;) {
                var u = Math.floor((o + s) / 2);
                n(t[u]) < i ? o = u + 1 : s = u
            }
            return o
        }, S.indexOf = u(1, S.findIndex, S.sortedIndex), S.lastIndexOf = u(-1, S.findLastIndex), S.range = function(t, e, n) {
            null == e && (e = t || 0, t = 0), n = n || 1;
            for (var r = Math.max(Math.ceil((e - t) / n), 0), i = Array(r), o = 0; o < r; o++, t += n) i[o] = t;
            return i
        };
        var N = function(t, e, n, r, i) {
            if (!(r instanceof e)) return t.apply(n, i);
            var o = E(t.prototype),
                s = t.apply(o, i);
            return S.isObject(s) ? s : o
        };
        S.bind = function(t, e) {
            if (b && t.bind === b) return b.apply(t, g.call(arguments, 1));
            if (!S.isFunction(t)) throw new TypeError("Bind must be called on a function");
            var n = g.call(arguments, 2);
            return function r() {
                return N(t, r, e, this, n.concat(g.call(arguments)))
            }
        }, S.partial = function(t) {
            var e = g.call(arguments, 1);
            return function n() {
                for (var r = 0, i = e.length, o = Array(i), s = 0; s < i; s++) o[s] = e[s] === S ? arguments[r++] : e[s];
                for (; r < arguments.length;) o.push(arguments[r++]);
                return N(t, n, this, this, o)
            }
        }, S.bindAll = function(t) {
            var e, n, r = arguments.length;
            if (r <= 1) throw new Error("bindAll must be passed function names");
            for (e = 1; e < r; e++) n = arguments[e], t[n] = S.bind(t[n], t);
            return t
        }, S.memoize = function(t, e) {
            var n = function r(n) {
                var i = r.cache,
                    o = "" + (e ? e.apply(this, arguments) : n);
                return S.has(i, o) || (i[o] = t.apply(this, arguments)), i[o]
            };
            return n.cache = {}, n
        }, S.delay = function(t, e) {
            var n = g.call(arguments, 2);
            return setTimeout(function() {
                return t.apply(null, n)
            }, e)
        }, S.defer = S.partial(S.delay, S, 1), S.throttle = function(t, e, n) {
            var r, i, o, s = null,
                u = 0;
            n || (n = {});
            var a = function() {
                u = n.leading === !1 ? 0 : S.now(), s = null, o = t.apply(r, i), s || (r = i = null)
            };
            return function() {
                var l = S.now();
                u || n.leading !== !1 || (u = l);
                var c = e - (l - u);
                return r = this, i = arguments, c <= 0 || c > e ? (s && (clearTimeout(s), s = null), u = l, o = t.apply(r, i), s || (r = i = null)) : s || n.trailing === !1 || (s = setTimeout(a, c)), o
            }
        }, S.debounce = function(t, e, n) {
            var r, i, o, s, u, a = function l() {
                var a = S.now() - s;
                a < e && a >= 0 ? r = setTimeout(l, e - a) : (r = null, n || (u = t.apply(o, i), r || (o = i = null)))
            };
            return function() {
                o = this, i = arguments, s = S.now();
                var l = n && !r;
                return r || (r = setTimeout(a, e)), l && (u = t.apply(o, i), o = i = null), u
            }
        }, S.wrap = function(t, e) {
            return S.partial(e, t)
        }, S.negate = function(t) {
            return function() {
                return !t.apply(this, arguments)
            }
        }, S.compose = function() {
            var t = arguments,
                e = t.length - 1;
            return function() {
                for (var n = e, r = t[e].apply(this, arguments); n--;) r = t[n].call(this, r);
                return r
            }
        }, S.after = function(t, e) {
            return function() {
                if (--t < 1) return e.apply(this, arguments)
            }
        }, S.before = function(t, e) {
            var n;
            return function() {
                return --t > 0 && (n = e.apply(this, arguments)), t <= 1 && (e = null), n
            }
        }, S.once = S.partial(S.before, 2);
        var L = !{
                toString: null
            }.propertyIsEnumerable("toString"),
            R = ["valueOf", "isPrototypeOf", "toString", "propertyIsEnumerable", "hasOwnProperty", "toLocaleString"];
        S.keys = function(t) {
            if (!S.isObject(t)) return [];
            if (_) return _(t);
            var e = [];
            for (var n in t) S.has(t, n) && e.push(n);
            return L && a(t, e), e
        }, S.allKeys = function(t) {
            if (!S.isObject(t)) return [];
            var e = [];
            for (var n in t) e.push(n);
            return L && a(t, e), e
        }, S.values = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = t[e[i]];
            return r
        }, S.mapObject = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = i.length, s = {}, u = 0; u < o; u++) r = i[u], s[r] = e(t[r], r, t);
            return s
        }, S.pairs = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = [e[i], t[e[i]]];
            return r
        }, S.invert = function(t) {
            for (var e = {}, n = S.keys(t), r = 0, i = n.length; r < i; r++) e[t[n[r]]] = n[r];
            return e
        }, S.functions = S.methods = function(t) {
            var e = [];
            for (var n in t) S.isFunction(t[n]) && e.push(n);
            return e.sort()
        }, S.extend = O(S.allKeys), S.extendOwn = S.assign = O(S.keys), S.findKey = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = 0, s = i.length; o < s; o++)
                if (r = i[o], e(t[r], r, t)) return r
        }, S.pick = function(t, e, n) {
            var r, i, o = {},
                s = t;
            if (null == s) return o;
            S.isFunction(e) ? (i = S.allKeys(s), r = k(e, n)) : (i = I(arguments, !1, !1, 1), r = function(t, e, n) {
                return e in n
            }, s = Object(s));
            for (var u = 0, a = i.length; u < a; u++) {
                var l = i[u],
                    c = s[l];
                r(c, l, s) && (o[l] = c)
            }
            return o
        }, S.omit = function(t, e, n) {
            if (S.isFunction(e)) e = S.negate(e);
            else {
                var r = S.map(I(arguments, !1, !1, 1), String);
                e = function(t, e) {
                    return !S.contains(r, e)
                }
            }
            return S.pick(t, e, n)
        }, S.defaults = O(S.allKeys, !0), S.create = function(t, e) {
            var n = E(t);
            return e && S.extendOwn(n, e), n
        }, S.clone = function(t) {
            return S.isObject(t) ? S.isArray(t) ? t.slice() : S.extend({}, t) : t
        }, S.tap = function(t, e) {
            return e(t), t
        }, S.isMatch = function(t, e) {
            var n = S.keys(e),
                r = n.length;
            if (null == t) return !r;
            for (var i = Object(t), o = 0; o < r; o++) {
                var s = n[o];
                if (e[s] !== i[s] || !(s in i)) return !1
            }
            return !0
        };
        var q = function X(t, e, n, r) {
            if (t === e) return 0 !== t || 1 / t === 1 / e;
            if (null == t || null == e) return t === e;
            t instanceof S && (t = t._wrapped), e instanceof S && (e = e._wrapped);
            var i = v.call(t);
            if (i !== v.call(e)) return !1;
            switch (i) {
                case "[object RegExp]":
                case "[object String]":
                    return "" + t == "" + e;
                case "[object Number]":
                    return +t !== +t ? +e !== +e : 0 === +t ? 1 / +t === 1 / e : +t === +e;
                case "[object Date]":
                case "[object Boolean]":
                    return +t === +e
            }
            var s = "[object Array]" === i;
            if (!s) {
                if ("object" != ("undefined" == typeof t ? "undefined" : o(t)) || "object" != ("undefined" == typeof e ? "undefined" : o(e))) return !1;
                var u = t.constructor,
                    a = e.constructor;
                if (u !== a && !(S.isFunction(u) && u instanceof u && S.isFunction(a) && a instanceof a) && "constructor" in t && "constructor" in e) return !1
            }
            n = n || [], r = r || [];
            for (var l = n.length; l--;)
                if (n[l] === t) return r[l] === e;
            if (n.push(t), r.push(e), s) {
                if (l = t.length, l !== e.length) return !1;
                for (; l--;)
                    if (!X(t[l], e[l], n, r)) return !1
            } else {
                var c, f = S.keys(t);
                if (l = f.length, S.keys(e).length !== l) return !1;
                for (; l--;)
                    if (c = f[l], !S.has(e, c) || !X(t[c], e[c], n, r)) return !1
            }
            return n.pop(), r.pop(), !0
        };
        S.isEqual = function(t, e) {
            return q(t, e)
        }, S.isEmpty = function(t) {
            return null == t || (C(t) && (S.isArray(t) || S.isString(t) || S.isArguments(t)) ? 0 === t.length : 0 === S.keys(t).length)
        }, S.isElement = function(t) {
            return !(!t || 1 !== t.nodeType)
        }, S.isArray = y || function(t) {
            return "[object Array]" === v.call(t)
        }, S.isObject = function(t) {
            var e = "undefined" == typeof t ? "undefined" : o(t);
            return "function" === e || "object" === e && !!t
        }, S.each(["Arguments", "Function", "String", "Number", "Date", "RegExp", "Error"], function(t) {
            S["is" + t] = function(e) {
                return v.call(e) === "[object " + t + "]"
            }
        }), S.isArguments(arguments) || (S.isArguments = function(t) {
            return S.has(t, "callee")
        }), "function" != typeof /./ && "object" != ("undefined" == typeof Int8Array ? "undefined" : o(Int8Array)) && (S.isFunction = function(t) {
            return "function" == typeof t || !1
        }), S.isFinite = function(t) {
            return isFinite(t) && !isNaN(parseFloat(t))
        }, S.isNaN = function(t) {
            return S.isNumber(t) && t !== +t
        }, S.isBoolean = function(t) {
            return t === !0 || t === !1 || "[object Boolean]" === v.call(t)
        }, S.isNull = function(t) {
            return null === t
        }, S.isUndefined = function(t) {
            return void 0 === t
        }, S.has = function(t, e) {
            return null != t && m.call(t, e)
        }, S.noConflict = function() {
            return l._ = c, this
        }, S.identity = function(t) {
            return t
        }, S.constant = function(t) {
            return function() {
                return t
            }
        }, S.noop = function() {}, S.property = M, S.propertyOf = function(t) {
            return null == t ? function() {} : function(e) {
                return t[e]
            }
        }, S.matcher = S.matches = function(t) {
            return t = S.extendOwn({}, t),
                function(e) {
                    return S.isMatch(e, t)
                }
        }, S.times = function(t, e, n) {
            var r = Array(Math.max(0, t));
            e = k(e, n, 1);
            for (var i = 0; i < t; i++) r[i] = e(i);
            return r
        }, S.random = function(t, e) {
            return null == e && (e = t, t = 0), t + Math.floor(Math.random() * (e - t + 1))
        }, S.now = Date.now || function() {
            return (new Date).getTime()
        };
        var F = {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': "&quot;",
                "'": "&#x27;",
                "`": "&#x60;"
            },
            P = S.invert(F),
            B = function(t) {
                var e = function(e) {
                        return t[e]
                    },
                    n = "(?:" + S.keys(t).join("|") + ")",
                    r = RegExp(n),
                    i = RegExp(n, "g");
                return function(t) {
                    return t = null == t ? "" : "" + t, r.test(t) ? t.replace(i, e) : t
                }
            };
        S.escape = B(F), S.unescape = B(P), S.result = function(t, e, n) {
            var r = null == t ? void 0 : t[e];
            return void 0 === r && (r = n), S.isFunction(r) ? r.call(t) : r
        };
        var W = 0;
        S.uniqueId = function(t) {
            var e = ++W + "";
            return t ? t + e : e
        }, S.templateSettings = {
            evaluate: /<%([\s\S]+?)%>/g,
            interpolate: /<%=([\s\S]+?)%>/g,
            escape: /<%-([\s\S]+?)%>/g
        };
        var D = /(.)^/,
            H = {
                "'": "'",
                "\\": "\\",
                "\r": "r",
                "\n": "n",
                "\u2028": "u2028",
                "\u2029": "u2029"
            },
            U = /\\|'|\r|\n|\u2028|\u2029/g,
            V = function(t) {
                return "\\" + H[t]
            };
        S.template = function(t, e, n) {
            !e && n && (e = n), e = S.defaults({}, e, S.templateSettings);
            var r = RegExp([(e.escape || D).source, (e.interpolate || D).source, (e.evaluate || D).source].join("|") + "|$", "g"),
                i = 0,
                o = "__p+='";
            t.replace(r, function(e, n, r, s, u) {
                return o += t.slice(i, u).replace(U, V), i = u + e.length, n ? o += "'+\n((__t=(" + n + "))==null?'':_.escape(__t))+\n'" : r ? o += "'+\n((__t=(" + r + "))==null?'':__t)+\n'" : s && (o += "';\n" + s + "\n__p+='"), e
            }), o += "';\n", e.variable || (o = "with(obj||{}){\n" + o + "}\n"), o = "var __t,__p='',__j=Array.prototype.join,print=function(){__p+=__j.call(arguments,'');};\n" + o + "return __p;\n";
            try {
                var s = new Function(e.variable || "obj", "_", o)
            } catch (u) {
                throw u.source = o, u
            }
            var a = function(t) {
                    return s.call(this, t, S)
                },
                l = e.variable || "obj";
            return a.source = "function(" + l + "){\n" + o + "}", a
        }, S.chain = function(t) {
            var e = S(t);
            return e._chain = !0, e
        };
        var $ = function(t, e) {
            return t._chain ? S(e).chain() : e
        };
        S.mixin = function(t) {
            S.each(S.functions(t), function(e) {
                var n = S[e] = t[e];
                S.prototype[e] = function() {
                    var t = [this._wrapped];
                    return p.apply(t, arguments), $(this, n.apply(S, t))
                }
            })
        }, S.mixin(S), S.each(["pop", "push", "reverse", "shift", "sort", "splice", "unshift"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                var n = this._wrapped;
                return e.apply(n, arguments), "shift" !== t && "splice" !== t || 0 !== n.length || delete n[0], $(this, n)
            }
        }), S.each(["concat", "join", "slice"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                return $(this, e.apply(this._wrapped, arguments))
            }
        }), S.prototype.value = function() {
            return this._wrapped
        }, S.prototype.valueOf = S.prototype.toJSON = S.prototype.value, S.prototype.toString = function() {
            return "" + this._wrapped
        }, r = [], i = function() {
            return S
        }.apply(e, r), !(void 0 !== i && (t.exports = i))
    }).call(void 0)
}, function(t, e, n) {
    "use strict";
    var r = n(22),
        i = n(23),
        o = n(42),
        s = n(41),
        u = [],
        a = u.slice,
        l = function(t, e) {
            e || (e = {}), e.model && (this.model = e.model), void 0 !== e.comparator && (this.comparator = e.comparator), this._reset(), this.initialize.apply(this, arguments), t && this.reset(t, o.extend({
                silent: !0
            }, e))
        },
        c = {
            add: !0,
            remove: !0,
            merge: !0
        },
        f = {
            add: !0,
            remove: !1
        };
    o.extend(l.prototype, r, {
        model: s,
        initialize: function() {},
        toJSON: function(t) {
            return this.map(function(e) {
                return e.toJSON(t)
            })
        },
        sync: function() {
            return Backbone.sync.apply(this, arguments)
        },
        add: function(t, e) {
            return this.set(t, o.extend({
                merge: !1
            }, e, f))
        },
        remove: function(t, e) {
            var n = !o.isArray(t);
            t = n ? [t] : o.clone(t), e || (e = {});
            for (var r = 0, i = t.length; r < i; r++) {
                var s = t[r] = this.get(t[r]);
                if (s) {
                    var u = this.modelId(s.attributes);
                    null != u && delete this._byId[u], delete this._byId[s.cid];
                    var a = this.indexOf(s);
                    this.models.splice(a, 1), this.length--, e.silent || (e.index = a, s.trigger("remove", s, this, e)), this._removeReference(s, e)
                }
            }
            return n ? t[0] : t
        },
        set: function(t, e) {
            e = o.defaults({}, e, c), e.parse && (t = this.parse(t, e));
            var n = !o.isArray(t);
            t = n ? t ? [t] : [] : t.slice();
            for (var r, i, s, u, a, l = e.at, f = this.comparator && null == l && e.sort !== !1, h = o.isString(this.comparator) ? this.comparator : null, d = [], p = [], g = {}, v = e.add, m = e.merge, y = e.remove, _ = !(f || !v || !y) && [], b = 0, x = t.length; b < x; b++) {
                if (s = t[b], u = this.get(s)) y && (g[u.cid] = !0), m && s !== u && (s = this._isModel(s) ? s.attributes : s, e.parse && (s = u.parse(s, e)), u.set(s, e), f && !a && u.hasChanged(h) && (a = !0)), t[b] = u;
                else if (v) {
                    if (i = t[b] = this._prepareModel(s, e), !i) continue;
                    d.push(i), this._addReference(i, e)
                }
                i = u || i, i && (r = this.modelId(i.attributes), !_ || !i.isNew() && g[r] || _.push(i), g[r] = !0)
            }
            if (y) {
                for (var b = 0, x = this.length; b < x; b++) g[(i = this.models[b]).cid] || p.push(i);
                p.length && this.remove(p, e)
            }
            if (d.length || _ && _.length)
                if (f && (a = !0), this.length += d.length, null != l)
                    for (var b = 0, x = d.length; b < x; b++) this.models.splice(l + b, 0, d[b]);
                else {
                    _ && (this.models.length = 0);
                    for (var w = _ || d, b = 0, x = w.length; b < x; b++) this.models.push(w[b])
                } if (a && this.sort({
                    silent: !0
                }), !e.silent) {
                for (var S = null != l ? o.clone(e) : e, b = 0, x = d.length; b < x; b++) null != l && (S.index = l + b), (i = d[b]).trigger("add", i, this, S);
                (a || _ && _.length) && this.trigger("sort", this, e)
            }
            return n ? t[0] : t
        },
        reset: function(t, e) {
            e || (e = {});
            for (var n = 0, r = this.models.length; n < r; n++) this._removeReference(this.models[n], e);
            return e.previousModels = this.models, this._reset(), t = this.add(t, o.extend({
                silent: !0
            }, e)), e.silent || this.trigger("reset", this, e), t
        },
        push: function(t, e) {
            return this.add(t, o.extend({
                at: this.length
            }, e))
        },
        pop: function(t) {
            var e = this.at(this.length - 1);
            return this.remove(e, t), e
        },
        unshift: function(t, e) {
            return this.add(t, o.extend({
                at: 0
            }, e))
        },
        shift: function(t) {
            var e = this.at(0);
            return this.remove(e, t), e
        },
        slice: function() {
            return a.apply(this.models, arguments)
        },
        get: function(t) {
            if (null != t) {
                var e = this.modelId(this._isModel(t) ? t.attributes : t);
                return this._byId[t] || this._byId[e] || this._byId[t.cid]
            }
        },
        at: function(t) {
            return t < 0 && (t += this.length), this.models[t]
        },
        where: function(t, e) {
            return o.isEmpty(t) ? e ? void 0 : [] : this[e ? "find" : "filter"](function(e) {
                for (var n in t)
                    if (t[n] !== e.get(n)) return !1;
                return !0
            })
        },
        findWhere: function(t) {
            return this.where(t, !0)
        },
        sort: function(t) {
            if (!this.comparator) throw new Error("Cannot sort a set without a comparator");
            return t || (t = {}), o.isString(this.comparator) || 1 === this.comparator.length ? this.models = this.sortBy(this.comparator, this) : this.models.sort(o.bind(this.comparator, this)), t.silent || this.trigger("sort", this, t), this
        },
        pluck: function(t) {
            return o.invoke(this.models, "get", t)
        },
        fetch: function(t) {
            t = t ? o.clone(t) : {}, void 0 === t.parse && (t.parse = !0);
            var e = t.success,
                n = this;
            return t.success = function(r) {
                n[t.reset ? "reset" : "set"](r, t), e && e(n, r, t), n.trigger("sync", n, r, t)
            }, wrapError(this, t), this.sync("read", this, t)
        },
        create: function(t, e) {
            if (e = e ? o.clone(e) : {}, !(t = this._prepareModel(t, e))) return !1;
            e.wait || this.add(t, e);
            var n = this,
                r = e.success;
            return e.success = function(t, i) {
                e.wait && n.add(t, e), r && r(t, i, e)
            }, t.save(null, e), t
        },
        parse: function(t, e) {
            return t
        },
        clone: function() {
            return new this.constructor(this.models, {
                model: this.model,
                comparator: this.comparator
            })
        },
        modelId: function(t) {
            return t[this.model.prototype.idAttribute || "id"]
        },
        _reset: function() {
            this.length = 0, this.models = [], this._byId = {}
        },
        _prepareModel: function(t, e) {
            if (this._isModel(t)) return t.collection || (t.collection = this), t;
            e = e ? o.clone(e) : {}, e.collection = this;
            var n = new this.model(t, e);
            return n.validationError ? (this.trigger("invalid", this, n.validationError, e), !1) : n
        },
        _isModel: function(t) {
            return t instanceof s
        },
        _addReference: function(t, e) {
            this._byId[t.cid] = t;
            var n = this.modelId(t.attributes);
            null != n && (this._byId[n] = t), t.on("all", this._onModelEvent, this)
        },
        _removeReference: function(t, e) {
            this === t.collection && delete t.collection, t.off("all", this._onModelEvent, this)
        },
        _onModelEvent: function(t, e, n, r) {
            if ("add" !== t && "remove" !== t || n === this) {
                if ("destroy" === t && this.remove(e, r), "change" === t) {
                    var i = this.modelId(e.previousAttributes()),
                        o = this.modelId(e.attributes);
                    i !== o && (null != i && delete this._byId[i], null != o && (this._byId[o] = e))
                }
                this.trigger.apply(this, arguments)
            }
        }
    });
    var h = ["forEach", "each", "map", "collect", "reduce", "foldl", "inject", "reduceRight", "foldr", "find", "detect", "filter", "select", "reject", "every", "all", "some", "any", "include", "contains", "invoke", "max", "min", "toArray", "size", "first", "head", "take", "initial", "rest", "tail", "drop", "last", "without", "difference", "indexOf", "shuffle", "lastIndexOf", "isEmpty", "chain", "sample", "partition"];
    o.each(h, function(t) {
        o[t] && (l.prototype[t] = function() {
            var e = a.call(arguments);
            return e.unshift(this.models), o[t].apply(o, e)
        })
    });
    var d = ["groupBy", "countBy", "sortBy", "indexBy"];
    o.each(d, function(t) {
        o[t] && (l.prototype[t] = function(e, n) {
            var r = o.isFunction(e) ? e : function(t) {
                return t.get(e)
            };
            return o[t](this.models, r, n)
        })
    }), l.extend = i, t.exports = l
}, function(t, e, n) {
    "use strict";
    var r = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    ! function() {
        function n() {
            return {
                keys: Object.keys || function(t) {
                    if ("object" !== ("undefined" == typeof t ? "undefined" : r(t)) && "function" != typeof t || null === t) throw new TypeError("keys() called on a non-object");
                    var e, n = [];
                    for (e in t) t.hasOwnProperty(e) && (n[n.length] = e);
                    return n
                },
                uniqueId: function(t) {
                    var e = ++a + "";
                    return t ? t + e : e
                },
                has: function(t, e) {
                    return s.call(t, e)
                },
                each: function(t, e, n) {
                    if (null != t)
                        if (o && t.forEach === o) t.forEach(e, n);
                        else if (t.length === +t.length)
                        for (var r = 0, i = t.length; r < i; r++) e.call(n, t[r], r, t);
                    else
                        for (var s in t) this.has(t, s) && e.call(n, t[s], s, t)
                },
                once: function(t) {
                    var e, n = !1;
                    return function() {
                        return n ? e : (n = !0, e = t.apply(this, arguments), t = null, e)
                    }
                }
            }
        }
        var i, o = Array.prototype.forEach,
            s = Object.prototype.hasOwnProperty,
            u = Array.prototype.slice,
            a = 0,
            l = n();
        i = {
            on: function(t, e, n) {
                return f(this, "on", t, [e, n]) && e ? (this._events || (this._events = {}), (this._events[t] || (this._events[t] = [])).push({
                    callback: e,
                    context: n,
                    ctx: n || this
                }), this) : this
            },
            once: function p(t, e, n) {
                if (!f(this, "once", t, [e, n]) || !e) return this;
                var r = this,
                    p = l.once(function() {
                        r.off(t, p), e.apply(this, arguments)
                    });
                return p._callback = e, this.on(t, p, n)
            },
            off: function(t, e, n) {
                var r, i, o, s, u, a, c, h;
                if (!this._events || !f(this, "off", t, [e, n])) return this;
                if (!t && !e && !n) return this._events = {}, this;
                for (s = t ? [t] : l.keys(this._events), u = 0, a = s.length; u < a; u++)
                    if (t = s[u], o = this._events[t]) {
                        if (this._events[t] = r = [], e || n)
                            for (c = 0, h = o.length; c < h; c++) i = o[c], (e && e !== i.callback && e !== i.callback._callback || n && n !== i.context) && r.push(i);
                        r.length || delete this._events[t]
                    } return this
            },
            trigger: function(t) {
                if (!this._events) return this;
                var e = u.call(arguments, 1);
                if (!f(this, "trigger", t, e)) return this;
                var n = this._events[t],
                    r = this._events.all;
                return n && h(n, e), r && h(r, arguments), this
            },
            stopListening: function(t, e, n) {
                var i = this._listeners;
                if (!i) return this;
                var o = !e && !n;
                "object" === ("undefined" == typeof e ? "undefined" : r(e)) && (n = this), t && ((i = {})[t._listenerId] = t);
                for (var s in i) i[s].off(e, n, this), o && delete this._listeners[s];
                return this
            }
        };
        var c = /\s+/,
            f = function(t, e, n, i) {
                if (!n) return !0;
                if ("object" === ("undefined" == typeof n ? "undefined" : r(n))) {
                    for (var o in n) t[e].apply(t, [o, n[o]].concat(i));
                    return !1
                }
                if (c.test(n)) {
                    for (var s = n.split(c), u = 0, a = s.length; u < a; u++) t[e].apply(t, [s[u]].concat(i));
                    return !1
                }
                return !0
            },
            h = function(t, e) {
                var n, r = -1,
                    i = t.length,
                    o = e[0],
                    s = e[1],
                    u = e[2];
                switch (e.length) {
                    case 0:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx);
                        return;
                    case 1:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o);
                        return;
                    case 2:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o, s);
                        return;
                    case 3:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o, s, u);
                        return;
                    default:
                        for (; ++r < i;)(n = t[r]).callback.apply(n.ctx, e)
                }
            },
            d = {
                listenTo: "on",
                listenToOnce: "once"
            };
        l.each(d, function(t, e) {
            i[e] = function(e, n, i) {
                return (this._listeners || (this._listeners = {}))[e._listenerId || (e._listenerId = l.uniqueId("l"))] = e, "object" === ("undefined" == typeof n ? "undefined" : r(n)) && (i = this), e[t](n, i, this), this
            }
        }), i.bind = i.on, i.unbind = i.off, i.mixin = function(t) {
            var e = ["on", "once", "off", "trigger", "stopListening", "listenTo", "listenToOnce", "bind", "unbind"];
            return l.each(e, function(e) {
                t[e] = this[e]
            }, this), t
        }, "undefined" != typeof t && t.exports && (e = t.exports = i), e.BackboneEvents = i
    }(void 0)
}, function(t, e, n) {
    "use strict";
    var r = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    ! function() {
        function n() {
            return {
                keys: Object.keys || function(t) {
                    if ("object" !== ("undefined" == typeof t ? "undefined" : r(t)) && "function" != typeof t || null === t) throw new TypeError("keys() called on a non-object");
                    var e, n = [];
                    for (e in t) t.hasOwnProperty(e) && (n[n.length] = e);
                    return n
                },
                uniqueId: function(t) {
                    var e = ++a + "";
                    return t ? t + e : e
                },
                has: function(t, e) {
                    return s.call(t, e)
                },
                each: function(t, e, n) {
                    if (null != t)
                        if (o && t.forEach === o) t.forEach(e, n);
                        else if (t.length === +t.length)
                        for (var r = 0, i = t.length; r < i; r++) e.call(n, t[r], r, t);
                    else
                        for (var s in t) this.has(t, s) && e.call(n, t[s], s, t)
                },
                once: function(t) {
                    var e, n = !1;
                    return function() {
                        return n ? e : (n = !0, e = t.apply(this, arguments), t = null, e)
                    }
                }
            }
        }
        var i, o = Array.prototype.forEach,
            s = Object.prototype.hasOwnProperty,
            u = Array.prototype.slice,
            a = 0,
            l = n();
        i = {
            on: function(t, e, n) {
                return f(this, "on", t, [e, n]) && e ? (this._events || (this._events = {}), (this._events[t] || (this._events[t] = [])).push({
                    callback: e,
                    context: n,
                    ctx: n || this
                }), this) : this
            },
            once: function p(t, e, n) {
                if (!f(this, "once", t, [e, n]) || !e) return this;
                var r = this,
                    p = l.once(function() {
                        r.off(t, p), e.apply(this, arguments)
                    });
                return p._callback = e, this.on(t, p, n)
            },
            off: function(t, e, n) {
                var r, i, o, s, u, a, c, h;
                if (!this._events || !f(this, "off", t, [e, n])) return this;
                if (!t && !e && !n) return this._events = {}, this;
                for (s = t ? [t] : l.keys(this._events), u = 0, a = s.length; u < a; u++)
                    if (t = s[u], o = this._events[t]) {
                        if (this._events[t] = r = [], e || n)
                            for (c = 0, h = o.length; c < h; c++) i = o[c], (e && e !== i.callback && e !== i.callback._callback || n && n !== i.context) && r.push(i);
                        r.length || delete this._events[t]
                    } return this
            },
            trigger: function(t) {
                if (!this._events) return this;
                var e = u.call(arguments, 1);
                if (!f(this, "trigger", t, e)) return this;
                var n = this._events[t],
                    r = this._events.all;
                return n && h(n, e), r && h(r, arguments), this
            },
            stopListening: function(t, e, n) {
                var i = this._listeners;
                if (!i) return this;
                var o = !e && !n;
                "object" === ("undefined" == typeof e ? "undefined" : r(e)) && (n = this), t && ((i = {})[t._listenerId] = t);
                for (var s in i) i[s].off(e, n, this), o && delete this._listeners[s];
                return this
            }
        };
        var c = /\s+/,
            f = function(t, e, n, i) {
                if (!n) return !0;
                if ("object" === ("undefined" == typeof n ? "undefined" : r(n))) {
                    for (var o in n) t[e].apply(t, [o, n[o]].concat(i));
                    return !1
                }
                if (c.test(n)) {
                    for (var s = n.split(c), u = 0, a = s.length; u < a; u++) t[e].apply(t, [s[u]].concat(i));
                    return !1
                }
                return !0
            },
            h = function(t, e) {
                var n, r = -1,
                    i = t.length,
                    o = e[0],
                    s = e[1],
                    u = e[2];
                switch (e.length) {
                    case 0:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx);
                        return;
                    case 1:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o);
                        return;
                    case 2:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o, s);
                        return;
                    case 3:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o, s, u);
                        return;
                    default:
                        for (; ++r < i;)(n = t[r]).callback.apply(n.ctx, e)
                }
            },
            d = {
                listenTo: "on",
                listenToOnce: "once"
            };
        l.each(d, function(t, e) {
            i[e] = function(e, n, i) {
                return (this._listeners || (this._listeners = {}))[e._listenerId || (e._listenerId = l.uniqueId("l"))] = e, "object" === ("undefined" == typeof n ? "undefined" : r(n)) && (i = this), e[t](n, i, this), this
            }
        }), i.bind = i.on, i.unbind = i.off, i.mixin = function(t) {
            var e = ["on", "once", "off", "trigger", "stopListening", "listenTo", "listenToOnce", "bind", "unbind"];
            return l.each(e, function(e) {
                t[e] = this[e]
            }, this), t
        }, "undefined" != typeof t && t.exports && (e = t.exports = i), e.BackboneEvents = i
    }(void 0)
}, function(t, e, n) {
    "use strict";
    t.exports = n(57)
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    ! function(s) {
        "object" === o(e) ? t.exports = s() : (r = s, i = "function" == typeof r ? r.call(e, n, e, t) : r, !(void 0 !== i && (t.exports = i)))
    }(function() {
        "use strict";
        var t = {
            has: function(t, e) {
                return Object.prototype.hasOwnProperty.call(t, e)
            },
            extend: function(t) {
                for (var e = 1; e < arguments.length; ++e) {
                    var n = arguments[e];
                    if (n)
                        for (var r in n) t[r] = n[r]
                }
                return t
            }
        };
        return function(e, n) {
            var r, i = this;
            r = e && t.has(e, "constructor") ? e.constructor : function() {
                return i.apply(this, arguments)
            }, t.extend(r, i, n);
            var o = function() {
                this.constructor = r
            };
            return o.prototype = i.prototype, r.prototype = new o, e && t.extend(r.prototype, e), r.__super__ = i.prototype, r
        }
    })
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    (function() {
        function n(t) {
            function e(e, n, r, i, o, s) {
                for (; o >= 0 && o < s; o += t) {
                    var u = i ? i[o] : o;
                    r = n(r, e[u], u, e)
                }
                return r
            }
            return function(n, r, i, o) {
                r = k(r, o, 4);
                var s = !C(n) && S.keys(n),
                    u = (s || n).length,
                    a = t > 0 ? 0 : u - 1;
                return arguments.length < 3 && (i = n[s ? s[a] : a], a += t), e(n, r, i, s, a, u)
            }
        }

        function s(t) {
            return function(e, n, r) {
                n = j(n, r);
                for (var i = A(e), o = t > 0 ? 0 : i - 1; o >= 0 && o < i; o += t)
                    if (n(e[o], o, e)) return o;
                return -1
            }
        }

        function u(t, e, n) {
            return function(r, i, o) {
                var s = 0,
                    u = A(r);
                if ("number" == typeof o) t > 0 ? s = o >= 0 ? o : Math.max(o + u, s) : u = o >= 0 ? Math.min(o + 1, u) : o + u + 1;
                else if (n && o && u) return o = n(r, i), r[o] === i ? o : -1;
                if (i !== i) return o = e(g.call(r, s, u), S.isNaN), o >= 0 ? o + s : -1;
                for (o = t > 0 ? s : u - 1; o >= 0 && o < u; o += t)
                    if (r[o] === i) return o;
                return -1
            }
        }

        function a(t, e) {
            var n = R.length,
                r = t.constructor,
                i = S.isFunction(r) && r.prototype || h,
                o = "constructor";
            for (S.has(t, o) && !S.contains(e, o) && e.push(o); n--;) o = R[n], o in t && t[o] !== i[o] && !S.contains(e, o) && e.push(o)
        }
        var l = this,
            c = l._,
            f = Array.prototype,
            h = Object.prototype,
            d = Function.prototype,
            p = f.push,
            g = f.slice,
            v = h.toString,
            m = h.hasOwnProperty,
            y = Array.isArray,
            _ = Object.keys,
            b = d.bind,
            x = Object.create,
            w = function() {},
            S = function G(t) {
                return t instanceof G ? t : this instanceof G ? void(this._wrapped = t) : new G(t)
            };
        "undefined" != typeof t && t.exports && (e = t.exports = S), e._ = S, S.VERSION = "1.8.3";
        var k = function(t, e, n) {
                if (void 0 === e) return t;
                switch (null == n ? 3 : n) {
                    case 1:
                        return function(n) {
                            return t.call(e, n)
                        };
                    case 2:
                        return function(n, r) {
                            return t.call(e, n, r)
                        };
                    case 3:
                        return function(n, r, i) {
                            return t.call(e, n, r, i)
                        };
                    case 4:
                        return function(n, r, i, o) {
                            return t.call(e, n, r, i, o)
                        }
                }
                return function() {
                    return t.apply(e, arguments)
                }
            },
            j = function(t, e, n) {
                return null == t ? S.identity : S.isFunction(t) ? k(t, e, n) : S.isObject(t) ? S.matcher(t) : S.property(t)
            };
        S.iteratee = function(t, e) {
            return j(t, e, 1 / 0)
        };
        var O = function(t, e) {
                return function(n) {
                    var r = arguments.length;
                    if (r < 2 || null == n) return n;
                    for (var i = 1; i < r; i++)
                        for (var o = arguments[i], s = t(o), u = s.length, a = 0; a < u; a++) {
                            var l = s[a];
                            e && void 0 !== n[l] || (n[l] = o[l])
                        }
                    return n
                }
            },
            E = function(t) {
                if (!S.isObject(t)) return {};
                if (x) return x(t);
                w.prototype = t;
                var e = new w;
                return w.prototype = null, e
            },
            M = function(t) {
                return function(e) {
                    return null == e ? void 0 : e[t]
                }
            },
            z = Math.pow(2, 53) - 1,
            A = M("length"),
            C = function(t) {
                var e = A(t);
                return "number" == typeof e && e >= 0 && e <= z
            };
        S.each = S.forEach = function(t, e, n) {
            e = k(e, n);
            var r, i;
            if (C(t))
                for (r = 0, i = t.length; r < i; r++) e(t[r], r, t);
            else {
                var o = S.keys(t);
                for (r = 0, i = o.length; r < i; r++) e(t[o[r]], o[r], t)
            }
            return t
        }, S.map = S.collect = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = Array(i), s = 0; s < i; s++) {
                var u = r ? r[s] : s;
                o[s] = e(t[u], u, t)
            }
            return o
        }, S.reduce = S.foldl = S.inject = n(1), S.reduceRight = S.foldr = n(-1), S.find = S.detect = function(t, e, n) {
            var r;
            if (r = C(t) ? S.findIndex(t, e, n) : S.findKey(t, e, n), void 0 !== r && r !== -1) return t[r]
        }, S.filter = S.select = function(t, e, n) {
            var r = [];
            return e = j(e, n), S.each(t, function(t, n, i) {
                e(t, n, i) && r.push(t)
            }), r
        }, S.reject = function(t, e, n) {
            return S.filter(t, S.negate(j(e)), n)
        }, S.every = S.all = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (!e(t[s], s, t)) return !1
            }
            return !0
        }, S.some = S.any = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (e(t[s], s, t)) return !0
            }
            return !1
        }, S.contains = S.includes = S.include = function(t, e, n, r) {
            return C(t) || (t = S.values(t)), ("number" != typeof n || r) && (n = 0), S.indexOf(t, e, n) >= 0
        }, S.invoke = function(t, e) {
            var n = g.call(arguments, 2),
                r = S.isFunction(e);
            return S.map(t, function(t) {
                var i = r ? e : t[e];
                return null == i ? i : i.apply(t, n)
            })
        }, S.pluck = function(t, e) {
            return S.map(t, S.property(e))
        }, S.where = function(t, e) {
            return S.filter(t, S.matcher(e))
        }, S.findWhere = function(t, e) {
            return S.find(t, S.matcher(e))
        }, S.max = function(t, e, n) {
            var r, i, o = -(1 / 0),
                s = -(1 / 0);
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r > o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i > s || i === -(1 / 0) && o === -(1 / 0)) && (o = t, s = i)
            });
            return o
        }, S.min = function(t, e, n) {
            var r, i, o = 1 / 0,
                s = 1 / 0;
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r < o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i < s || i === 1 / 0 && o === 1 / 0) && (o = t, s = i)
            });
            return o
        }, S.shuffle = function(t) {
            for (var e, n = C(t) ? t : S.values(t), r = n.length, i = Array(r), o = 0; o < r; o++) e = S.random(0, o), e !== o && (i[o] = i[e]), i[e] = n[o];
            return i
        }, S.sample = function(t, e, n) {
            return null == e || n ? (C(t) || (t = S.values(t)), t[S.random(t.length - 1)]) : S.shuffle(t).slice(0, Math.max(0, e))
        }, S.sortBy = function(t, e, n) {
            return e = j(e, n), S.pluck(S.map(t, function(t, n, r) {
                return {
                    value: t,
                    index: n,
                    criteria: e(t, n, r)
                }
            }).sort(function(t, e) {
                var n = t.criteria,
                    r = e.criteria;
                if (n !== r) {
                    if (n > r || void 0 === n) return 1;
                    if (n < r || void 0 === r) return -1
                }
                return t.index - e.index
            }), "value")
        };
        var T = function(t) {
            return function(e, n, r) {
                var i = {};
                return n = j(n, r), S.each(e, function(r, o) {
                    t(i, r, n(r, o, e))
                }), i
            }
        };
        S.groupBy = T(function(t, e, n) {
            S.has(t, n) ? t[n].push(e) : t[n] = [e]
        }), S.indexBy = T(function(t, e, n) {
            t[n] = e
        }), S.countBy = T(function(t, e, n) {
            S.has(t, n) ? t[n]++ : t[n] = 1
        }), S.toArray = function(t) {
            return t ? S.isArray(t) ? g.call(t) : C(t) ? S.map(t, S.identity) : S.values(t) : []
        }, S.size = function(t) {
            return null == t ? 0 : C(t) ? t.length : S.keys(t).length
        }, S.partition = function(t, e, n) {
            e = j(e, n);
            var r = [],
                i = [];
            return S.each(t, function(t, n, o) {
                (e(t, n, o) ? r : i).push(t)
            }), [r, i]
        }, S.first = S.head = S.take = function(t, e, n) {
            if (null != t) return null == e || n ? t[0] : S.initial(t, t.length - e)
        }, S.initial = function(t, e, n) {
            return g.call(t, 0, Math.max(0, t.length - (null == e || n ? 1 : e)))
        }, S.last = function(t, e, n) {
            if (null != t) return null == e || n ? t[t.length - 1] : S.rest(t, Math.max(0, t.length - e))
        }, S.rest = S.tail = S.drop = function(t, e, n) {
            return g.call(t, null == e || n ? 1 : e)
        }, S.compact = function(t) {
            return S.filter(t, S.identity)
        };
        var I = function K(t, e, n, r) {
            for (var i = [], o = 0, s = r || 0, u = A(t); s < u; s++) {
                var a = t[s];
                if (C(a) && (S.isArray(a) || S.isArguments(a))) {
                    e || (a = K(a, e, n));
                    var l = 0,
                        c = a.length;
                    for (i.length += c; l < c;) i[o++] = a[l++]
                } else n || (i[o++] = a)
            }
            return i
        };
        S.flatten = function(t, e) {
            return I(t, e, !1)
        }, S.without = function(t) {
            return S.difference(t, g.call(arguments, 1))
        }, S.uniq = S.unique = function(t, e, n, r) {
            S.isBoolean(e) || (r = n, n = e, e = !1), null != n && (n = j(n, r));
            for (var i = [], o = [], s = 0, u = A(t); s < u; s++) {
                var a = t[s],
                    l = n ? n(a, s, t) : a;
                e ? (s && o === l || i.push(a), o = l) : n ? S.contains(o, l) || (o.push(l), i.push(a)) : S.contains(i, a) || i.push(a)
            }
            return i
        }, S.union = function() {
            return S.uniq(I(arguments, !0, !0))
        }, S.intersection = function(t) {
            for (var e = [], n = arguments.length, r = 0, i = A(t); r < i; r++) {
                var o = t[r];
                if (!S.contains(e, o)) {
                    for (var s = 1; s < n && S.contains(arguments[s], o); s++);
                    s === n && e.push(o)
                }
            }
            return e
        }, S.difference = function(t) {
            var e = I(arguments, !0, !0, 1);
            return S.filter(t, function(t) {
                return !S.contains(e, t)
            })
        }, S.zip = function() {
            return S.unzip(arguments)
        }, S.unzip = function(t) {
            for (var e = t && S.max(t, A).length || 0, n = Array(e), r = 0; r < e; r++) n[r] = S.pluck(t, r);
            return n
        }, S.object = function(t, e) {
            for (var n = {}, r = 0, i = A(t); r < i; r++) e ? n[t[r]] = e[r] : n[t[r][0]] = t[r][1];
            return n
        }, S.findIndex = s(1), S.findLastIndex = s(-1), S.sortedIndex = function(t, e, n, r) {
            n = j(n, r, 1);
            for (var i = n(e), o = 0, s = A(t); o < s;) {
                var u = Math.floor((o + s) / 2);
                n(t[u]) < i ? o = u + 1 : s = u
            }
            return o
        }, S.indexOf = u(1, S.findIndex, S.sortedIndex), S.lastIndexOf = u(-1, S.findLastIndex), S.range = function(t, e, n) {
            null == e && (e = t || 0, t = 0), n = n || 1;
            for (var r = Math.max(Math.ceil((e - t) / n), 0), i = Array(r), o = 0; o < r; o++, t += n) i[o] = t;
            return i
        };
        var N = function(t, e, n, r, i) {
            if (!(r instanceof e)) return t.apply(n, i);
            var o = E(t.prototype),
                s = t.apply(o, i);
            return S.isObject(s) ? s : o
        };
        S.bind = function(t, e) {
            if (b && t.bind === b) return b.apply(t, g.call(arguments, 1));
            if (!S.isFunction(t)) throw new TypeError("Bind must be called on a function");
            var n = g.call(arguments, 2);
            return function r() {
                return N(t, r, e, this, n.concat(g.call(arguments)))
            }
        }, S.partial = function(t) {
            var e = g.call(arguments, 1);
            return function n() {
                for (var r = 0, i = e.length, o = Array(i), s = 0; s < i; s++) o[s] = e[s] === S ? arguments[r++] : e[s];
                for (; r < arguments.length;) o.push(arguments[r++]);
                return N(t, n, this, this, o)
            }
        }, S.bindAll = function(t) {
            var e, n, r = arguments.length;
            if (r <= 1) throw new Error("bindAll must be passed function names");
            for (e = 1; e < r; e++) n = arguments[e], t[n] = S.bind(t[n], t);
            return t
        }, S.memoize = function(t, e) {
            var n = function r(n) {
                var i = r.cache,
                    o = "" + (e ? e.apply(this, arguments) : n);
                return S.has(i, o) || (i[o] = t.apply(this, arguments)), i[o]
            };
            return n.cache = {}, n
        }, S.delay = function(t, e) {
            var n = g.call(arguments, 2);
            return setTimeout(function() {
                return t.apply(null, n)
            }, e)
        }, S.defer = S.partial(S.delay, S, 1), S.throttle = function(t, e, n) {
            var r, i, o, s = null,
                u = 0;
            n || (n = {});
            var a = function() {
                u = n.leading === !1 ? 0 : S.now(), s = null, o = t.apply(r, i), s || (r = i = null)
            };
            return function() {
                var l = S.now();
                u || n.leading !== !1 || (u = l);
                var c = e - (l - u);
                return r = this, i = arguments, c <= 0 || c > e ? (s && (clearTimeout(s), s = null), u = l, o = t.apply(r, i), s || (r = i = null)) : s || n.trailing === !1 || (s = setTimeout(a, c)), o
            }
        }, S.debounce = function(t, e, n) {
            var r, i, o, s, u, a = function l() {
                var a = S.now() - s;
                a < e && a >= 0 ? r = setTimeout(l, e - a) : (r = null, n || (u = t.apply(o, i), r || (o = i = null)))
            };
            return function() {
                o = this, i = arguments, s = S.now();
                var l = n && !r;
                return r || (r = setTimeout(a, e)), l && (u = t.apply(o, i), o = i = null), u
            }
        }, S.wrap = function(t, e) {
            return S.partial(e, t)
        }, S.negate = function(t) {
            return function() {
                return !t.apply(this, arguments)
            }
        }, S.compose = function() {
            var t = arguments,
                e = t.length - 1;
            return function() {
                for (var n = e, r = t[e].apply(this, arguments); n--;) r = t[n].call(this, r);
                return r
            }
        }, S.after = function(t, e) {
            return function() {
                if (--t < 1) return e.apply(this, arguments)
            }
        }, S.before = function(t, e) {
            var n;
            return function() {
                return --t > 0 && (n = e.apply(this, arguments)), t <= 1 && (e = null), n
            }
        }, S.once = S.partial(S.before, 2);
        var L = !{
                toString: null
            }.propertyIsEnumerable("toString"),
            R = ["valueOf", "isPrototypeOf", "toString", "propertyIsEnumerable", "hasOwnProperty", "toLocaleString"];
        S.keys = function(t) {
            if (!S.isObject(t)) return [];
            if (_) return _(t);
            var e = [];
            for (var n in t) S.has(t, n) && e.push(n);
            return L && a(t, e), e
        }, S.allKeys = function(t) {
            if (!S.isObject(t)) return [];
            var e = [];
            for (var n in t) e.push(n);
            return L && a(t, e), e
        }, S.values = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = t[e[i]];
            return r
        }, S.mapObject = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = i.length, s = {}, u = 0; u < o; u++) r = i[u], s[r] = e(t[r], r, t);
            return s
        }, S.pairs = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = [e[i], t[e[i]]];
            return r
        }, S.invert = function(t) {
            for (var e = {}, n = S.keys(t), r = 0, i = n.length; r < i; r++) e[t[n[r]]] = n[r];
            return e
        }, S.functions = S.methods = function(t) {
            var e = [];
            for (var n in t) S.isFunction(t[n]) && e.push(n);
            return e.sort()
        }, S.extend = O(S.allKeys), S.extendOwn = S.assign = O(S.keys), S.findKey = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = 0, s = i.length; o < s; o++)
                if (r = i[o], e(t[r], r, t)) return r
        }, S.pick = function(t, e, n) {
            var r, i, o = {},
                s = t;
            if (null == s) return o;
            S.isFunction(e) ? (i = S.allKeys(s), r = k(e, n)) : (i = I(arguments, !1, !1, 1), r = function(t, e, n) {
                return e in n
            }, s = Object(s));
            for (var u = 0, a = i.length; u < a; u++) {
                var l = i[u],
                    c = s[l];
                r(c, l, s) && (o[l] = c)
            }
            return o
        }, S.omit = function(t, e, n) {
            if (S.isFunction(e)) e = S.negate(e);
            else {
                var r = S.map(I(arguments, !1, !1, 1), String);
                e = function(t, e) {
                    return !S.contains(r, e)
                }
            }
            return S.pick(t, e, n)
        }, S.defaults = O(S.allKeys, !0), S.create = function(t, e) {
            var n = E(t);
            return e && S.extendOwn(n, e), n
        }, S.clone = function(t) {
            return S.isObject(t) ? S.isArray(t) ? t.slice() : S.extend({}, t) : t
        }, S.tap = function(t, e) {
            return e(t), t
        }, S.isMatch = function(t, e) {
            var n = S.keys(e),
                r = n.length;
            if (null == t) return !r;
            for (var i = Object(t), o = 0; o < r; o++) {
                var s = n[o];
                if (e[s] !== i[s] || !(s in i)) return !1
            }
            return !0
        };
        var q = function X(t, e, n, r) {
            if (t === e) return 0 !== t || 1 / t === 1 / e;
            if (null == t || null == e) return t === e;
            t instanceof S && (t = t._wrapped), e instanceof S && (e = e._wrapped);
            var i = v.call(t);
            if (i !== v.call(e)) return !1;
            switch (i) {
                case "[object RegExp]":
                case "[object String]":
                    return "" + t == "" + e;
                case "[object Number]":
                    return +t !== +t ? +e !== +e : 0 === +t ? 1 / +t === 1 / e : +t === +e;
                case "[object Date]":
                case "[object Boolean]":
                    return +t === +e
            }
            var s = "[object Array]" === i;
            if (!s) {
                if ("object" != ("undefined" == typeof t ? "undefined" : o(t)) || "object" != ("undefined" == typeof e ? "undefined" : o(e))) return !1;
                var u = t.constructor,
                    a = e.constructor;
                if (u !== a && !(S.isFunction(u) && u instanceof u && S.isFunction(a) && a instanceof a) && "constructor" in t && "constructor" in e) return !1
            }
            n = n || [], r = r || [];
            for (var l = n.length; l--;)
                if (n[l] === t) return r[l] === e;
            if (n.push(t), r.push(e), s) {
                if (l = t.length, l !== e.length) return !1;
                for (; l--;)
                    if (!X(t[l], e[l], n, r)) return !1
            } else {
                var c, f = S.keys(t);
                if (l = f.length, S.keys(e).length !== l) return !1;
                for (; l--;)
                    if (c = f[l], !S.has(e, c) || !X(t[c], e[c], n, r)) return !1
            }
            return n.pop(), r.pop(), !0
        };
        S.isEqual = function(t, e) {
            return q(t, e)
        }, S.isEmpty = function(t) {
            return null == t || (C(t) && (S.isArray(t) || S.isString(t) || S.isArguments(t)) ? 0 === t.length : 0 === S.keys(t).length)
        }, S.isElement = function(t) {
            return !(!t || 1 !== t.nodeType)
        }, S.isArray = y || function(t) {
            return "[object Array]" === v.call(t)
        }, S.isObject = function(t) {
            var e = "undefined" == typeof t ? "undefined" : o(t);
            return "function" === e || "object" === e && !!t
        }, S.each(["Arguments", "Function", "String", "Number", "Date", "RegExp", "Error"], function(t) {
            S["is" + t] = function(e) {
                return v.call(e) === "[object " + t + "]"
            }
        }), S.isArguments(arguments) || (S.isArguments = function(t) {
            return S.has(t, "callee")
        }), "function" != typeof /./ && "object" != ("undefined" == typeof Int8Array ? "undefined" : o(Int8Array)) && (S.isFunction = function(t) {
            return "function" == typeof t || !1
        }), S.isFinite = function(t) {
            return isFinite(t) && !isNaN(parseFloat(t))
        }, S.isNaN = function(t) {
            return S.isNumber(t) && t !== +t
        }, S.isBoolean = function(t) {
            return t === !0 || t === !1 || "[object Boolean]" === v.call(t)
        }, S.isNull = function(t) {
            return null === t
        }, S.isUndefined = function(t) {
            return void 0 === t
        }, S.has = function(t, e) {
            return null != t && m.call(t, e)
        }, S.noConflict = function() {
            return l._ = c, this
        }, S.identity = function(t) {
            return t
        }, S.constant = function(t) {
            return function() {
                return t
            }
        }, S.noop = function() {}, S.property = M, S.propertyOf = function(t) {
            return null == t ? function() {} : function(e) {
                return t[e]
            }
        }, S.matcher = S.matches = function(t) {
            return t = S.extendOwn({}, t),
                function(e) {
                    return S.isMatch(e, t)
                }
        }, S.times = function(t, e, n) {
            var r = Array(Math.max(0, t));
            e = k(e, n, 1);
            for (var i = 0; i < t; i++) r[i] = e(i);
            return r
        }, S.random = function(t, e) {
            return null == e && (e = t, t = 0), t + Math.floor(Math.random() * (e - t + 1))
        }, S.now = Date.now || function() {
            return (new Date).getTime()
        };
        var F = {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': "&quot;",
                "'": "&#x27;",
                "`": "&#x60;"
            },
            P = S.invert(F),
            B = function(t) {
                var e = function(e) {
                        return t[e]
                    },
                    n = "(?:" + S.keys(t).join("|") + ")",
                    r = RegExp(n),
                    i = RegExp(n, "g");
                return function(t) {
                    return t = null == t ? "" : "" + t, r.test(t) ? t.replace(i, e) : t
                }
            };
        S.escape = B(F), S.unescape = B(P), S.result = function(t, e, n) {
            var r = null == t ? void 0 : t[e];
            return void 0 === r && (r = n), S.isFunction(r) ? r.call(t) : r
        };
        var W = 0;
        S.uniqueId = function(t) {
            var e = ++W + "";
            return t ? t + e : e
        }, S.templateSettings = {
            evaluate: /<%([\s\S]+?)%>/g,
            interpolate: /<%=([\s\S]+?)%>/g,
            escape: /<%-([\s\S]+?)%>/g
        };
        var D = /(.)^/,
            H = {
                "'": "'",
                "\\": "\\",
                "\r": "r",
                "\n": "n",
                "\u2028": "u2028",
                "\u2029": "u2029"
            },
            U = /\\|'|\r|\n|\u2028|\u2029/g,
            V = function(t) {
                return "\\" + H[t]
            };
        S.template = function(t, e, n) {
            !e && n && (e = n), e = S.defaults({}, e, S.templateSettings);
            var r = RegExp([(e.escape || D).source, (e.interpolate || D).source, (e.evaluate || D).source].join("|") + "|$", "g"),
                i = 0,
                o = "__p+='";
            t.replace(r, function(e, n, r, s, u) {
                return o += t.slice(i, u).replace(U, V), i = u + e.length, n ? o += "'+\n((__t=(" + n + "))==null?'':_.escape(__t))+\n'" : r ? o += "'+\n((__t=(" + r + "))==null?'':__t)+\n'" : s && (o += "';\n" + s + "\n__p+='"), e
            }), o += "';\n", e.variable || (o = "with(obj||{}){\n" + o + "}\n"), o = "var __t,__p='',__j=Array.prototype.join,print=function(){__p+=__j.call(arguments,'');};\n" + o + "return __p;\n";
            try {
                var s = new Function(e.variable || "obj", "_", o)
            } catch (u) {
                throw u.source = o, u
            }
            var a = function(t) {
                    return s.call(this, t, S)
                },
                l = e.variable || "obj";
            return a.source = "function(" + l + "){\n" + o + "}", a
        }, S.chain = function(t) {
            var e = S(t);
            return e._chain = !0, e
        };
        var $ = function(t, e) {
            return t._chain ? S(e).chain() : e
        };
        S.mixin = function(t) {
            S.each(S.functions(t), function(e) {
                var n = S[e] = t[e];
                S.prototype[e] = function() {
                    var t = [this._wrapped];
                    return p.apply(t, arguments), $(this, n.apply(S, t))
                }
            })
        }, S.mixin(S), S.each(["pop", "push", "reverse", "shift", "sort", "splice", "unshift"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                var n = this._wrapped;
                return e.apply(n, arguments), "shift" !== t && "splice" !== t || 0 !== n.length || delete n[0], $(this, n)
            }
        }), S.each(["concat", "join", "slice"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                return $(this, e.apply(this._wrapped, arguments))
            }
        }), S.prototype.value = function() {
            return this._wrapped
        }, S.prototype.valueOf = S.prototype.toJSON = S.prototype.value, S.prototype.toString = function() {
            return "" + this._wrapped
        }, r = [], i = function() {
            return S
        }.apply(e, r), !(void 0 !== i && (t.exports = i))
    }).call(void 0)
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(11),
        o = r(i),
        s = n(24),
        u = r(s),
        a = void 0;
    e["default"] = a = {
        parse: function(t) {
            var e = [];
            if ("[object Array]" === Object.prototype.toString.call(t)) var n = t;
            else var n = t.split("\n");
            if (n[0].slice(0, 6) === !1) throw new Error("Invalid CLUSTAL Header");
            for (var r = 0, i = 1, o = 0; r < n.length;) {
                r++;
                var s = n[r];
                if (null != s && 0 !== s.length)
                    if (0 !== s.trim().length) {
                        if (!u["default"].contains(s, "*")) {
                            1 === i && (o = 0, i = 0);
                            var a = /^(?:\s*)(\S+)(?:\s+)(\S+)(?:\s*)(\d*)(?:\s*|$)/g,
                                l = a.exec(s);
                            if (null != l) {
                                var c = l[1].trim(),
                                    f = l[2].trim();
                                if (o >= e.length) {
                                    var h = u["default"].getMeta(c.trim());
                                    c = h.name;
                                    var d = new u["default"].model(f, c, o);
                                    d.ids = h.ids || {}, d.details = h.details || {};
                                    var p = Object.keys(d.ids);
                                    p.length > 0 && (d.id = d.ids[p[0]]), e.push(d)
                                } else e[o].seq += f;
                                o++
                            } else console.log("parse error", s)
                        }
                    } else i = 1;
                else i = 1
            }
            return e
        }
    }, o["default"].mixin(a)
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(24),
        o = r(i),
        s = n(69),
        u = r(s),
        a = n(11),
        l = r(a),
        c = void 0;
    e["default"] = c = {
        getMeta: o["default"].getMeta,
        extend: function(t) {
            var e = (0, u["default"])(c);
            return c.getMeta = t, e
        },
        parse: function(t) {
            var e = [];
            if (!t || 0 === t.length) return [];
            "[object Array]" !== Object.prototype.toString.call(t) && (t = t.split("\n"));
            for (var n = c, r = n.getMeta, i = 0; i < t.length; i++) {
                var s = t[i];
                if (">" === s[0] || ";" === s[0]) {
                    var u = s.slice(1).trim(),
                        a = r(u.trim());
                    u = a.name;
                    var l = a.id || e.length,
                        f = new o["default"].model("", a.name, l);
                    f.ids = a.ids || {}, f.details = a.details || {}, e.push(f)
                } else f.seq += s
            }
            return e
        },
        write: function(t, e) {
            for (var n = "", r = 0; r < t.length; r++) {
                var i = t[r];
                null != e && (i = e(i)), n += ">" + i.name + "\n", n += o["default"].splitNChars(i.seq, 80).join("\n"), n += "\n"
            }
            return n
        }
    }, l["default"].mixin(c)
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(11),
        o = r(i),
        s = n(43),
        u = r(s),
        a = n(64),
        l = r(a),
        c = function() {};
    o["default"].mixin(c), e["default"] = c, c.parseLines = function(t) {
        var e = t.split("\n"),
            n = {},
            r = [];
        n.type = c._guessType(e);
        var i = 0;
        if ("jalview" === n.type) {
            var o = l["default"].readHeader(e);
            i = o.offset, n.colors = o.colors, r = o.features
        }
        for (var s = i; s < e.length; s++) {
            var u = e[s];
            0 !== u.length && "#" !== u[0] && (u = c.parseLine(u), void 0 !== u && r.push(u))
        }
        return {
            features: r,
            config: n
        }
    }, c._guessType = function(t) {
        return "##gff-version 3" === t[0].substring(0, 15) ? "gff3" : t[0].indexOf("#") < 0 && 2 === t[0].split("\t").length ? "jalview" : "gff3"
    }, c.parseSeqs = c.parse = function(t) {
        var e = c.parseLines(t),
            n = {};
        return e.features.forEach(function(t) {
            var e = t.seqname;
            void 0 === n[e] && (n[e] = []), delete t.seqname, n[e].push(t)
        }), delete e.features, e.seqs = n, e
    }, c.parseLine = function(t) {
        var e = {},
            n = t.split(/\s+/);
        if (1 !== n.length) {
            e.seqname = n[0], e.source = n[1], e.feature = n[2], e.start = parseInt(n[3]), e.end = parseInt(n[4]), e.score = n[5], e.strand = n[6], e.frame = n[7];
            var r = n.slice(8).join(" ");
            return Object.keys(e).forEach(function(t) {
                "string" == typeof e[t] && (e[t] = e[t].trim()), "." === e[t] && (e[t] = void 0)
            }), e.score && (e.score = parseFloat(e.score)), e.frame && (e.frame = parseInt(e.frame)), e.attributes = u["default"].extractKeys(r), e
        }
    }, c.exportLine = function(t) {
        var e = Object.keys(t.attributes).map(function(e) {
                return e + "=" + t.attributes[e]
            }).join(";"),
            n = [t.seqname, t.source, t.feature, t.start, t.end, t.score, t.strand, t.frame, e];
        return n = n.map(function(t) {
            return void 0 === t ? "." : t
        }), n.join("\t")
    }, c.exportLines = function(t) {
        return "##gff-version 3\n" + t.map(c.exportLine).join("\n")
    }, c.exportSeqs = c["export"] = function(t) {
        var e = [],
            n = function(t) {
                t.seqname = r, e.push(t)
            };
        for (var r in t) t[r].forEach(n);
        return c.exportLines(e)
    }, c.utils = u["default"]
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(43),
        i = {};
    e["default"] = i, i.readHeader = function(t) {
        for (var e, n = {}, r = 0, o = []; r < t.length; r++) {
            var s = t[r];
            if (s.indexOf("#") >= 0) break;
            var u = s.split(/\t/),
                a = u[0].trim();
            if ("GFF" === a) break;
            if (2 === u.length)
                if ("startgroup" === a) e = u[1].trim();
                else {
                    if ("endgroup" === a) {
                        e = "";
                        continue
                    }
                    n[u[0]] = i.parseColor(u[1])
                }
            else if (u.length >= 5) {
                var l = i.parseLine(u);
                e && (l.attributes.Parent = e), o.push(l)
            }
        }
        return {
            offset: r,
            colors: n,
            features: o
        }
    }, i.parseColor = function(t) {
        return t.indexOf(",") >= 0 ? (0, r.rgbToHex)(t.split(",").map(function(t) {
            return parseInt(t)
        })) : 6 === t.length && parseInt(t.charAt(0), 16) <= 16 && "bisque" !== t ? "#" + t : t
    }, i.parseLine = function(t) {
        var e = {
            attributes: {}
        };
        return e.attributes.Name = t[0].trim(), e.seqname = t[1].trim(), e.start = parseInt(t[3]), e.end = parseInt(t[4]), e.feature = t[5].trim(), "ID_NOT_SPECIFIED" === e.seqname && (e.seqname = t[2].trim()), e
    }
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    var i = n(11),
        o = r(i),
        s = function u(t) {
            return this.constructor != u ? new u(t) : (this.matrix = {}, this.parsingOrder = [], void 0 != t && this.parse(t), this)
        };
    o["default"].mixin(s), t.exports = s, s.prototype.parse = function(t) {
        return t.split("\n").forEach(function(t) {
            this.parseLine(t)
        }.bind(this)), this.buildMatrix(), this.matrix
    }, s.read = function(t, e) {
        return (new s).read(t, e)
    }, s.parse = function(t) {
        return (new s).parse(t)
    }, s.prototype.parseLine = function(t) {
        var e = t.charAt(0);
        if ("#" !== e) {
            this.parsingOrder.push(e);
            for (var n = t.substring(1), r = n.split(/\s+/).filter(function(t) {
                    return t.length > 0
                }).map(function(t) {
                    return parseInt(t)
                }), i = {}, o = 0; o < r.length; o++) i[this.parsingOrder[o]] = r[o];
            this.matrix[e] = i
        }
    }, s.prototype["export"] = function() {
        return s["export"](this.matrix)
    }, s["export"] = function(t) {
        var e = [],
            n = 1;
        "matrix" in t && (t = t.matrix);
        for (var r in t) {
            for (var i = r, o = Object.keys(t[r]), s = 0; s < n; s++) i += "\t" + t[r][o[s]];
            e.push(i), n++
        }
        return e.join("\n")
    }, s.prototype.buildMatrix = function() {
        for (var t = this.parsingOrder.length - 1; t >= 0; t--) {
            var e = this.parsingOrder[t],
                n = this.matrix[e];
            for (var r in n) this.matrix[r][e] = n[r]
        }
    }
}, function(t, e) {
    "use strict";

    function n(t) {
        for (var e = [], n = {}, i = t.split(/\s*(;|\(|\)|\[|\]|,|:|=)\s*/), o = 0; o < i.length; o++) {
            var s = i[o],
                u = void 0;
            switch (s) {
                case "(":
                    u = {}, n.children = [u], e.push(n), n = u;
                    break;
                case ",":
                    u = {}, e[e.length - 1].children.push(u), n = u;
                    break;
                case ")":
                    n = e.pop();
                    break;
                case ":":
                    break;
                default:
                    var a = i[o - 1];
                    if (")" == a || "(" == a || "," == a) n.name = s;
                    else if (":" == a) "undefined" == typeof s ? "undefined" : r(s), isNaN(s) || (n.branch_length = parseFloat(s));
                    else if ("=" == a) {
                        var l = i[o - 2];
                        switch (l) {
                            case "D":
                                n.duplication = s;
                                break;
                            case "G":
                                n.gene_id = s;
                                break;
                            case "T":
                                n.taxon_id = s
                        }
                    }
            }
        }
        return n
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    e["default"] = n
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e.parseNhx = e.parse = void 0;
    var i = n(68),
        o = r(i),
        s = n(66),
        u = r(s),
        a = {};
    a.parse = o["default"], a.parseNhx = u["default"], e["default"] = a, e.parse = o["default"], e.parseNhx = u["default"]
}, function(t, e) {
    "use strict";

    function n(t) {
        for (var e = [], n = {}, r = t.split(/\s*(;|\(|\)|,|:)\s*/), i = 0; i < r.length; i++) {
            var o = r[i],
                s = void 0;
            switch (o) {
                case "(":
                    s = {}, n.children = [s], e.push(n), n = s;
                    break;
                case ",":
                    s = {}, e[e.length - 1].children.push(s), n = s;
                    break;
                case ")":
                    n = e.pop();
                    break;
                case ":":
                    break;
                default:
                    var u = r[i - 1];
                    ")" == u || "(" == u || "," == u ? n.name = o : ":" == u && (n.branch_length = parseFloat(o))
            }
        }
        return n
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e["default"] = n
}, function(t, e) {
    "use strict";

    function n(t, e, n) {
        for (var r = [], i = t < e, o = n ? i ? e + 1 : e - 1 : e, s = t; i ? s < o : s > o; i ? s++ : s--) r.push(s);
        return r
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e["default"] = function(t) {
        t = t || {};
        for (var e = n(0, arguments.length, !1), r = 0; r < e.length; r++) {
            var i = e[r];
            if (arguments[i])
                for (var o = 0; o < arguments[i].length; o++) {
                    var s = arguments[i][o];
                    arguments[i].hasOwnProperty(s) && (t[s] = arguments[i][s])
                }
        }
        return t
    }
}, function(t, e, n) {
    "use strict";
    var r = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    ! function() {
        function n() {
            return {
                keys: Object.keys || function(t) {
                    if ("object" !== ("undefined" == typeof t ? "undefined" : r(t)) && "function" != typeof t || null === t) throw new TypeError("keys() called on a non-object");
                    var e, n = [];
                    for (e in t) t.hasOwnProperty(e) && (n[n.length] = e);
                    return n
                },
                uniqueId: function(t) {
                    var e = ++a + "";
                    return t ? t + e : e
                },
                has: function(t, e) {
                    return s.call(t, e)
                },
                each: function(t, e, n) {
                    if (null != t)
                        if (o && t.forEach === o) t.forEach(e, n);
                        else if (t.length === +t.length)
                        for (var r = 0, i = t.length; r < i; r++) e.call(n, t[r], r, t);
                    else
                        for (var s in t) this.has(t, s) && e.call(n, t[s], s, t)
                },
                once: function(t) {
                    var e, n = !1;
                    return function() {
                        return n ? e : (n = !0, e = t.apply(this, arguments), t = null, e)
                    }
                }
            }
        }
        var i, o = Array.prototype.forEach,
            s = Object.prototype.hasOwnProperty,
            u = Array.prototype.slice,
            a = 0,
            l = n();
        i = {
            on: function(t, e, n) {
                return f(this, "on", t, [e, n]) && e ? (this._events || (this._events = {}), (this._events[t] || (this._events[t] = [])).push({
                    callback: e,
                    context: n,
                    ctx: n || this
                }), this) : this
            },
            once: function p(t, e, n) {
                if (!f(this, "once", t, [e, n]) || !e) return this;
                var r = this,
                    p = l.once(function() {
                        r.off(t, p), e.apply(this, arguments)
                    });
                return p._callback = e, this.on(t, p, n)
            },
            off: function(t, e, n) {
                var r, i, o, s, u, a, c, h;
                if (!this._events || !f(this, "off", t, [e, n])) return this;
                if (!t && !e && !n) return this._events = {}, this;
                for (s = t ? [t] : l.keys(this._events), u = 0, a = s.length; u < a; u++)
                    if (t = s[u], o = this._events[t]) {
                        if (this._events[t] = r = [], e || n)
                            for (c = 0, h = o.length; c < h; c++) i = o[c], (e && e !== i.callback && e !== i.callback._callback || n && n !== i.context) && r.push(i);
                        r.length || delete this._events[t]
                    } return this
            },
            trigger: function(t) {
                if (!this._events) return this;
                var e = u.call(arguments, 1);
                if (!f(this, "trigger", t, e)) return this;
                var n = this._events[t],
                    r = this._events.all;
                return n && h(n, e), r && h(r, arguments), this
            },
            stopListening: function(t, e, n) {
                var i = this._listeners;
                if (!i) return this;
                var o = !e && !n;
                "object" === ("undefined" == typeof e ? "undefined" : r(e)) && (n = this), t && ((i = {})[t._listenerId] = t);
                for (var s in i) i[s].off(e, n, this), o && delete this._listeners[s];
                return this
            }
        };
        var c = /\s+/,
            f = function(t, e, n, i) {
                if (!n) return !0;
                if ("object" === ("undefined" == typeof n ? "undefined" : r(n))) {
                    for (var o in n) t[e].apply(t, [o, n[o]].concat(i));
                    return !1
                }
                if (c.test(n)) {
                    for (var s = n.split(c), u = 0, a = s.length; u < a; u++) t[e].apply(t, [s[u]].concat(i));
                    return !1
                }
                return !0
            },
            h = function(t, e) {
                var n, r = -1,
                    i = t.length,
                    o = e[0],
                    s = e[1],
                    u = e[2];
                switch (e.length) {
                    case 0:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx);
                        return;
                    case 1:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o);
                        return;
                    case 2:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o, s);
                        return;
                    case 3:
                        for (; ++r < i;)(n = t[r]).callback.call(n.ctx, o, s, u);
                        return;
                    default:
                        for (; ++r < i;)(n = t[r]).callback.apply(n.ctx, e)
                }
            },
            d = {
                listenTo: "on",
                listenToOnce: "once"
            };
        l.each(d, function(t, e) {
            i[e] = function(e, n, i) {
                return (this._listeners || (this._listeners = {}))[e._listenerId || (e._listenerId = l.uniqueId("l"))] = e, "object" === ("undefined" == typeof n ? "undefined" : r(n)) && (i = this), e[t](n, i, this), this
            }
        }), i.bind = i.on, i.unbind = i.off, i.mixin = function(t) {
            var e = ["on", "once", "off", "trigger", "stopListening", "listenTo", "listenToOnce", "bind", "unbind"];
            return l.each(e, function(e) {
                t[e] = this[e]
            }, this), t
        }, "undefined" != typeof t && t.exports && (e = t.exports = i), e.BackboneEvents = i
    }(void 0)
}, function(t, e, n) {
    "use strict";
    t.exports = n(70)
}, function(t, e, n) {
    "use strict";
    t.exports.seq = n(73)
}, function(t, e) {
    "use strict";
    t.exports = function(t, e, n) {
        this.seq = t, this.name = e, this.id = n, this.meta = {}
    }
}, function(t, e, n) {
    "use strict";
    t.exports = n(81)
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    (function() {
        function n(t) {
            function e(e, n, r, i, o, s) {
                for (; o >= 0 && o < s; o += t) {
                    var u = i ? i[o] : o;
                    r = n(r, e[u], u, e)
                }
                return r
            }
            return function(n, r, i, o) {
                r = k(r, o, 4);
                var s = !C(n) && S.keys(n),
                    u = (s || n).length,
                    a = t > 0 ? 0 : u - 1;
                return arguments.length < 3 && (i = n[s ? s[a] : a], a += t), e(n, r, i, s, a, u)
            }
        }

        function s(t) {
            return function(e, n, r) {
                n = j(n, r);
                for (var i = A(e), o = t > 0 ? 0 : i - 1; o >= 0 && o < i; o += t)
                    if (n(e[o], o, e)) return o;
                return -1
            }
        }

        function u(t, e, n) {
            return function(r, i, o) {
                var s = 0,
                    u = A(r);
                if ("number" == typeof o) t > 0 ? s = o >= 0 ? o : Math.max(o + u, s) : u = o >= 0 ? Math.min(o + 1, u) : o + u + 1;
                else if (n && o && u) return o = n(r, i), r[o] === i ? o : -1;
                if (i !== i) return o = e(g.call(r, s, u), S.isNaN), o >= 0 ? o + s : -1;
                for (o = t > 0 ? s : u - 1; o >= 0 && o < u; o += t)
                    if (r[o] === i) return o;
                return -1
            }
        }

        function a(t, e) {
            var n = R.length,
                r = t.constructor,
                i = S.isFunction(r) && r.prototype || h,
                o = "constructor";
            for (S.has(t, o) && !S.contains(e, o) && e.push(o); n--;) o = R[n], o in t && t[o] !== i[o] && !S.contains(e, o) && e.push(o)
        }
        var l = this,
            c = l._,
            f = Array.prototype,
            h = Object.prototype,
            d = Function.prototype,
            p = f.push,
            g = f.slice,
            v = h.toString,
            m = h.hasOwnProperty,
            y = Array.isArray,
            _ = Object.keys,
            b = d.bind,
            x = Object.create,
            w = function() {},
            S = function G(t) {
                return t instanceof G ? t : this instanceof G ? void(this._wrapped = t) : new G(t)
            };
        "undefined" != typeof t && t.exports && (e = t.exports = S), e._ = S, S.VERSION = "1.8.3";
        var k = function(t, e, n) {
                if (void 0 === e) return t;
                switch (null == n ? 3 : n) {
                    case 1:
                        return function(n) {
                            return t.call(e, n)
                        };
                    case 2:
                        return function(n, r) {
                            return t.call(e, n, r)
                        };
                    case 3:
                        return function(n, r, i) {
                            return t.call(e, n, r, i)
                        };
                    case 4:
                        return function(n, r, i, o) {
                            return t.call(e, n, r, i, o)
                        }
                }
                return function() {
                    return t.apply(e, arguments)
                }
            },
            j = function(t, e, n) {
                return null == t ? S.identity : S.isFunction(t) ? k(t, e, n) : S.isObject(t) ? S.matcher(t) : S.property(t)
            };
        S.iteratee = function(t, e) {
            return j(t, e, 1 / 0)
        };
        var O = function(t, e) {
                return function(n) {
                    var r = arguments.length;
                    if (r < 2 || null == n) return n;
                    for (var i = 1; i < r; i++)
                        for (var o = arguments[i], s = t(o), u = s.length, a = 0; a < u; a++) {
                            var l = s[a];
                            e && void 0 !== n[l] || (n[l] = o[l])
                        }
                    return n
                }
            },
            E = function(t) {
                if (!S.isObject(t)) return {};
                if (x) return x(t);
                w.prototype = t;
                var e = new w;
                return w.prototype = null, e
            },
            M = function(t) {
                return function(e) {
                    return null == e ? void 0 : e[t]
                }
            },
            z = Math.pow(2, 53) - 1,
            A = M("length"),
            C = function(t) {
                var e = A(t);
                return "number" == typeof e && e >= 0 && e <= z
            };
        S.each = S.forEach = function(t, e, n) {
            e = k(e, n);
            var r, i;
            if (C(t))
                for (r = 0, i = t.length; r < i; r++) e(t[r], r, t);
            else {
                var o = S.keys(t);
                for (r = 0, i = o.length; r < i; r++) e(t[o[r]], o[r], t)
            }
            return t
        }, S.map = S.collect = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = Array(i), s = 0; s < i; s++) {
                var u = r ? r[s] : s;
                o[s] = e(t[u], u, t)
            }
            return o
        }, S.reduce = S.foldl = S.inject = n(1), S.reduceRight = S.foldr = n(-1), S.find = S.detect = function(t, e, n) {
            var r;
            if (r = C(t) ? S.findIndex(t, e, n) : S.findKey(t, e, n), void 0 !== r && r !== -1) return t[r]
        }, S.filter = S.select = function(t, e, n) {
            var r = [];
            return e = j(e, n), S.each(t, function(t, n, i) {
                e(t, n, i) && r.push(t)
            }), r
        }, S.reject = function(t, e, n) {
            return S.filter(t, S.negate(j(e)), n)
        }, S.every = S.all = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (!e(t[s], s, t)) return !1
            }
            return !0
        }, S.some = S.any = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (e(t[s], s, t)) return !0
            }
            return !1
        }, S.contains = S.includes = S.include = function(t, e, n, r) {
            return C(t) || (t = S.values(t)), ("number" != typeof n || r) && (n = 0), S.indexOf(t, e, n) >= 0
        }, S.invoke = function(t, e) {
            var n = g.call(arguments, 2),
                r = S.isFunction(e);
            return S.map(t, function(t) {
                var i = r ? e : t[e];
                return null == i ? i : i.apply(t, n)
            })
        }, S.pluck = function(t, e) {
            return S.map(t, S.property(e))
        }, S.where = function(t, e) {
            return S.filter(t, S.matcher(e))
        }, S.findWhere = function(t, e) {
            return S.find(t, S.matcher(e))
        }, S.max = function(t, e, n) {
            var r, i, o = -(1 / 0),
                s = -(1 / 0);
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r > o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i > s || i === -(1 / 0) && o === -(1 / 0)) && (o = t, s = i)
            });
            return o
        }, S.min = function(t, e, n) {
            var r, i, o = 1 / 0,
                s = 1 / 0;
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r < o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i < s || i === 1 / 0 && o === 1 / 0) && (o = t, s = i)
            });
            return o
        }, S.shuffle = function(t) {
            for (var e, n = C(t) ? t : S.values(t), r = n.length, i = Array(r), o = 0; o < r; o++) e = S.random(0, o), e !== o && (i[o] = i[e]), i[e] = n[o];
            return i
        }, S.sample = function(t, e, n) {
            return null == e || n ? (C(t) || (t = S.values(t)), t[S.random(t.length - 1)]) : S.shuffle(t).slice(0, Math.max(0, e))
        }, S.sortBy = function(t, e, n) {
            return e = j(e, n), S.pluck(S.map(t, function(t, n, r) {
                return {
                    value: t,
                    index: n,
                    criteria: e(t, n, r)
                }
            }).sort(function(t, e) {
                var n = t.criteria,
                    r = e.criteria;
                if (n !== r) {
                    if (n > r || void 0 === n) return 1;
                    if (n < r || void 0 === r) return -1
                }
                return t.index - e.index
            }), "value")
        };
        var T = function(t) {
            return function(e, n, r) {
                var i = {};
                return n = j(n, r), S.each(e, function(r, o) {
                    t(i, r, n(r, o, e))
                }), i
            }
        };
        S.groupBy = T(function(t, e, n) {
            S.has(t, n) ? t[n].push(e) : t[n] = [e]
        }), S.indexBy = T(function(t, e, n) {
            t[n] = e
        }), S.countBy = T(function(t, e, n) {
            S.has(t, n) ? t[n]++ : t[n] = 1
        }), S.toArray = function(t) {
            return t ? S.isArray(t) ? g.call(t) : C(t) ? S.map(t, S.identity) : S.values(t) : []
        }, S.size = function(t) {
            return null == t ? 0 : C(t) ? t.length : S.keys(t).length
        }, S.partition = function(t, e, n) {
            e = j(e, n);
            var r = [],
                i = [];
            return S.each(t, function(t, n, o) {
                (e(t, n, o) ? r : i).push(t)
            }), [r, i]
        }, S.first = S.head = S.take = function(t, e, n) {
            if (null != t) return null == e || n ? t[0] : S.initial(t, t.length - e)
        }, S.initial = function(t, e, n) {
            return g.call(t, 0, Math.max(0, t.length - (null == e || n ? 1 : e)))
        }, S.last = function(t, e, n) {
            if (null != t) return null == e || n ? t[t.length - 1] : S.rest(t, Math.max(0, t.length - e))
        }, S.rest = S.tail = S.drop = function(t, e, n) {
            return g.call(t, null == e || n ? 1 : e)
        }, S.compact = function(t) {
            return S.filter(t, S.identity)
        };
        var I = function K(t, e, n, r) {
            for (var i = [], o = 0, s = r || 0, u = A(t); s < u; s++) {
                var a = t[s];
                if (C(a) && (S.isArray(a) || S.isArguments(a))) {
                    e || (a = K(a, e, n));
                    var l = 0,
                        c = a.length;
                    for (i.length += c; l < c;) i[o++] = a[l++]
                } else n || (i[o++] = a)
            }
            return i
        };
        S.flatten = function(t, e) {
            return I(t, e, !1)
        }, S.without = function(t) {
            return S.difference(t, g.call(arguments, 1))
        }, S.uniq = S.unique = function(t, e, n, r) {
            S.isBoolean(e) || (r = n, n = e, e = !1), null != n && (n = j(n, r));
            for (var i = [], o = [], s = 0, u = A(t); s < u; s++) {
                var a = t[s],
                    l = n ? n(a, s, t) : a;
                e ? (s && o === l || i.push(a), o = l) : n ? S.contains(o, l) || (o.push(l), i.push(a)) : S.contains(i, a) || i.push(a)
            }
            return i
        }, S.union = function() {
            return S.uniq(I(arguments, !0, !0))
        }, S.intersection = function(t) {
            for (var e = [], n = arguments.length, r = 0, i = A(t); r < i; r++) {
                var o = t[r];
                if (!S.contains(e, o)) {
                    for (var s = 1; s < n && S.contains(arguments[s], o); s++);
                    s === n && e.push(o)
                }
            }
            return e
        }, S.difference = function(t) {
            var e = I(arguments, !0, !0, 1);
            return S.filter(t, function(t) {
                return !S.contains(e, t)
            })
        }, S.zip = function() {
            return S.unzip(arguments)
        }, S.unzip = function(t) {
            for (var e = t && S.max(t, A).length || 0, n = Array(e), r = 0; r < e; r++) n[r] = S.pluck(t, r);
            return n
        }, S.object = function(t, e) {
            for (var n = {}, r = 0, i = A(t); r < i; r++) e ? n[t[r]] = e[r] : n[t[r][0]] = t[r][1];
            return n
        }, S.findIndex = s(1), S.findLastIndex = s(-1), S.sortedIndex = function(t, e, n, r) {
            n = j(n, r, 1);
            for (var i = n(e), o = 0, s = A(t); o < s;) {
                var u = Math.floor((o + s) / 2);
                n(t[u]) < i ? o = u + 1 : s = u
            }
            return o
        }, S.indexOf = u(1, S.findIndex, S.sortedIndex), S.lastIndexOf = u(-1, S.findLastIndex), S.range = function(t, e, n) {
            null == e && (e = t || 0, t = 0), n = n || 1;
            for (var r = Math.max(Math.ceil((e - t) / n), 0), i = Array(r), o = 0; o < r; o++, t += n) i[o] = t;
            return i
        };
        var N = function(t, e, n, r, i) {
            if (!(r instanceof e)) return t.apply(n, i);
            var o = E(t.prototype),
                s = t.apply(o, i);
            return S.isObject(s) ? s : o
        };
        S.bind = function(t, e) {
            if (b && t.bind === b) return b.apply(t, g.call(arguments, 1));
            if (!S.isFunction(t)) throw new TypeError("Bind must be called on a function");
            var n = g.call(arguments, 2);
            return function r() {
                return N(t, r, e, this, n.concat(g.call(arguments)))
            }
        }, S.partial = function(t) {
            var e = g.call(arguments, 1);
            return function n() {
                for (var r = 0, i = e.length, o = Array(i), s = 0; s < i; s++) o[s] = e[s] === S ? arguments[r++] : e[s];
                for (; r < arguments.length;) o.push(arguments[r++]);
                return N(t, n, this, this, o)
            }
        }, S.bindAll = function(t) {
            var e, n, r = arguments.length;
            if (r <= 1) throw new Error("bindAll must be passed function names");
            for (e = 1; e < r; e++) n = arguments[e], t[n] = S.bind(t[n], t);
            return t
        }, S.memoize = function(t, e) {
            var n = function r(n) {
                var i = r.cache,
                    o = "" + (e ? e.apply(this, arguments) : n);
                return S.has(i, o) || (i[o] = t.apply(this, arguments)), i[o]
            };
            return n.cache = {}, n
        }, S.delay = function(t, e) {
            var n = g.call(arguments, 2);
            return setTimeout(function() {
                return t.apply(null, n)
            }, e)
        }, S.defer = S.partial(S.delay, S, 1), S.throttle = function(t, e, n) {
            var r, i, o, s = null,
                u = 0;
            n || (n = {});
            var a = function() {
                u = n.leading === !1 ? 0 : S.now(), s = null, o = t.apply(r, i), s || (r = i = null)
            };
            return function() {
                var l = S.now();
                u || n.leading !== !1 || (u = l);
                var c = e - (l - u);
                return r = this, i = arguments, c <= 0 || c > e ? (s && (clearTimeout(s), s = null), u = l, o = t.apply(r, i), s || (r = i = null)) : s || n.trailing === !1 || (s = setTimeout(a, c)), o
            }
        }, S.debounce = function(t, e, n) {
            var r, i, o, s, u, a = function l() {
                var a = S.now() - s;
                a < e && a >= 0 ? r = setTimeout(l, e - a) : (r = null, n || (u = t.apply(o, i), r || (o = i = null)))
            };
            return function() {
                o = this, i = arguments, s = S.now();
                var l = n && !r;
                return r || (r = setTimeout(a, e)), l && (u = t.apply(o, i), o = i = null), u
            }
        }, S.wrap = function(t, e) {
            return S.partial(e, t)
        }, S.negate = function(t) {
            return function() {
                return !t.apply(this, arguments)
            }
        }, S.compose = function() {
            var t = arguments,
                e = t.length - 1;
            return function() {
                for (var n = e, r = t[e].apply(this, arguments); n--;) r = t[n].call(this, r);
                return r
            }
        }, S.after = function(t, e) {
            return function() {
                if (--t < 1) return e.apply(this, arguments)
            }
        }, S.before = function(t, e) {
            var n;
            return function() {
                return --t > 0 && (n = e.apply(this, arguments)), t <= 1 && (e = null), n
            }
        }, S.once = S.partial(S.before, 2);
        var L = !{
                toString: null
            }.propertyIsEnumerable("toString"),
            R = ["valueOf", "isPrototypeOf", "toString", "propertyIsEnumerable", "hasOwnProperty", "toLocaleString"];
        S.keys = function(t) {
            if (!S.isObject(t)) return [];
            if (_) return _(t);
            var e = [];
            for (var n in t) S.has(t, n) && e.push(n);
            return L && a(t, e), e
        }, S.allKeys = function(t) {
            if (!S.isObject(t)) return [];
            var e = [];
            for (var n in t) e.push(n);
            return L && a(t, e), e
        }, S.values = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = t[e[i]];
            return r
        }, S.mapObject = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = i.length, s = {}, u = 0; u < o; u++) r = i[u], s[r] = e(t[r], r, t);
            return s
        }, S.pairs = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = [e[i], t[e[i]]];
            return r
        }, S.invert = function(t) {
            for (var e = {}, n = S.keys(t), r = 0, i = n.length; r < i; r++) e[t[n[r]]] = n[r];
            return e
        }, S.functions = S.methods = function(t) {
            var e = [];
            for (var n in t) S.isFunction(t[n]) && e.push(n);
            return e.sort()
        }, S.extend = O(S.allKeys), S.extendOwn = S.assign = O(S.keys), S.findKey = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = 0, s = i.length; o < s; o++)
                if (r = i[o], e(t[r], r, t)) return r
        }, S.pick = function(t, e, n) {
            var r, i, o = {},
                s = t;
            if (null == s) return o;
            S.isFunction(e) ? (i = S.allKeys(s), r = k(e, n)) : (i = I(arguments, !1, !1, 1), r = function(t, e, n) {
                return e in n
            }, s = Object(s));
            for (var u = 0, a = i.length; u < a; u++) {
                var l = i[u],
                    c = s[l];
                r(c, l, s) && (o[l] = c)
            }
            return o
        }, S.omit = function(t, e, n) {
            if (S.isFunction(e)) e = S.negate(e);
            else {
                var r = S.map(I(arguments, !1, !1, 1), String);
                e = function(t, e) {
                    return !S.contains(r, e)
                }
            }
            return S.pick(t, e, n)
        }, S.defaults = O(S.allKeys, !0), S.create = function(t, e) {
            var n = E(t);
            return e && S.extendOwn(n, e), n
        }, S.clone = function(t) {
            return S.isObject(t) ? S.isArray(t) ? t.slice() : S.extend({}, t) : t
        }, S.tap = function(t, e) {
            return e(t), t
        }, S.isMatch = function(t, e) {
            var n = S.keys(e),
                r = n.length;
            if (null == t) return !r;
            for (var i = Object(t), o = 0; o < r; o++) {
                var s = n[o];
                if (e[s] !== i[s] || !(s in i)) return !1
            }
            return !0
        };
        var q = function X(t, e, n, r) {
            if (t === e) return 0 !== t || 1 / t === 1 / e;
            if (null == t || null == e) return t === e;
            t instanceof S && (t = t._wrapped), e instanceof S && (e = e._wrapped);
            var i = v.call(t);
            if (i !== v.call(e)) return !1;
            switch (i) {
                case "[object RegExp]":
                case "[object String]":
                    return "" + t == "" + e;
                case "[object Number]":
                    return +t !== +t ? +e !== +e : 0 === +t ? 1 / +t === 1 / e : +t === +e;
                case "[object Date]":
                case "[object Boolean]":
                    return +t === +e
            }
            var s = "[object Array]" === i;
            if (!s) {
                if ("object" != ("undefined" == typeof t ? "undefined" : o(t)) || "object" != ("undefined" == typeof e ? "undefined" : o(e))) return !1;
                var u = t.constructor,
                    a = e.constructor;
                if (u !== a && !(S.isFunction(u) && u instanceof u && S.isFunction(a) && a instanceof a) && "constructor" in t && "constructor" in e) return !1
            }
            n = n || [], r = r || [];
            for (var l = n.length; l--;)
                if (n[l] === t) return r[l] === e;
            if (n.push(t), r.push(e), s) {
                if (l = t.length, l !== e.length) return !1;
                for (; l--;)
                    if (!X(t[l], e[l], n, r)) return !1
            } else {
                var c, f = S.keys(t);
                if (l = f.length, S.keys(e).length !== l) return !1;
                for (; l--;)
                    if (c = f[l], !S.has(e, c) || !X(t[c], e[c], n, r)) return !1
            }
            return n.pop(), r.pop(), !0
        };
        S.isEqual = function(t, e) {
            return q(t, e)
        }, S.isEmpty = function(t) {
            return null == t || (C(t) && (S.isArray(t) || S.isString(t) || S.isArguments(t)) ? 0 === t.length : 0 === S.keys(t).length)
        }, S.isElement = function(t) {
            return !(!t || 1 !== t.nodeType)
        }, S.isArray = y || function(t) {
            return "[object Array]" === v.call(t)
        }, S.isObject = function(t) {
            var e = "undefined" == typeof t ? "undefined" : o(t);
            return "function" === e || "object" === e && !!t
        }, S.each(["Arguments", "Function", "String", "Number", "Date", "RegExp", "Error"], function(t) {
            S["is" + t] = function(e) {
                return v.call(e) === "[object " + t + "]"
            }
        }), S.isArguments(arguments) || (S.isArguments = function(t) {
            return S.has(t, "callee")
        }), "function" != typeof /./ && "object" != ("undefined" == typeof Int8Array ? "undefined" : o(Int8Array)) && (S.isFunction = function(t) {
            return "function" == typeof t || !1
        }), S.isFinite = function(t) {
            return isFinite(t) && !isNaN(parseFloat(t))
        }, S.isNaN = function(t) {
            return S.isNumber(t) && t !== +t
        }, S.isBoolean = function(t) {
            return t === !0 || t === !1 || "[object Boolean]" === v.call(t)
        }, S.isNull = function(t) {
            return null === t
        }, S.isUndefined = function(t) {
            return void 0 === t
        }, S.has = function(t, e) {
            return null != t && m.call(t, e)
        }, S.noConflict = function() {
            return l._ = c, this
        }, S.identity = function(t) {
            return t
        }, S.constant = function(t) {
            return function() {
                return t
            }
        }, S.noop = function() {}, S.property = M, S.propertyOf = function(t) {
            return null == t ? function() {} : function(e) {
                return t[e]
            }
        }, S.matcher = S.matches = function(t) {
            return t = S.extendOwn({}, t),
                function(e) {
                    return S.isMatch(e, t)
                }
        }, S.times = function(t, e, n) {
            var r = Array(Math.max(0, t));
            e = k(e, n, 1);
            for (var i = 0; i < t; i++) r[i] = e(i);
            return r
        }, S.random = function(t, e) {
            return null == e && (e = t, t = 0), t + Math.floor(Math.random() * (e - t + 1))
        }, S.now = Date.now || function() {
            return (new Date).getTime()
        };
        var F = {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': "&quot;",
                "'": "&#x27;",
                "`": "&#x60;"
            },
            P = S.invert(F),
            B = function(t) {
                var e = function(e) {
                        return t[e]
                    },
                    n = "(?:" + S.keys(t).join("|") + ")",
                    r = RegExp(n),
                    i = RegExp(n, "g");
                return function(t) {
                    return t = null == t ? "" : "" + t, r.test(t) ? t.replace(i, e) : t
                }
            };
        S.escape = B(F), S.unescape = B(P), S.result = function(t, e, n) {
            var r = null == t ? void 0 : t[e];
            return void 0 === r && (r = n), S.isFunction(r) ? r.call(t) : r
        };
        var W = 0;
        S.uniqueId = function(t) {
            var e = ++W + "";
            return t ? t + e : e
        }, S.templateSettings = {
            evaluate: /<%([\s\S]+?)%>/g,
            interpolate: /<%=([\s\S]+?)%>/g,
            escape: /<%-([\s\S]+?)%>/g
        };
        var D = /(.)^/,
            H = {
                "'": "'",
                "\\": "\\",
                "\r": "r",
                "\n": "n",
                "\u2028": "u2028",
                "\u2029": "u2029"
            },
            U = /\\|'|\r|\n|\u2028|\u2029/g,
            V = function(t) {
                return "\\" + H[t]
            };
        S.template = function(t, e, n) {
            !e && n && (e = n), e = S.defaults({}, e, S.templateSettings);
            var r = RegExp([(e.escape || D).source, (e.interpolate || D).source, (e.evaluate || D).source].join("|") + "|$", "g"),
                i = 0,
                o = "__p+='";
            t.replace(r, function(e, n, r, s, u) {
                return o += t.slice(i, u).replace(U, V), i = u + e.length, n ? o += "'+\n((__t=(" + n + "))==null?'':_.escape(__t))+\n'" : r ? o += "'+\n((__t=(" + r + "))==null?'':__t)+\n'" : s && (o += "';\n" + s + "\n__p+='"), e
            }), o += "';\n", e.variable || (o = "with(obj||{}){\n" + o + "}\n"), o = "var __t,__p='',__j=Array.prototype.join,print=function(){__p+=__j.call(arguments,'');};\n" + o + "return __p;\n";
            try {
                var s = new Function(e.variable || "obj", "_", o)
            } catch (u) {
                throw u.source = o, u
            }
            var a = function(t) {
                    return s.call(this, t, S)
                },
                l = e.variable || "obj";
            return a.source = "function(" + l + "){\n" + o + "}", a
        }, S.chain = function(t) {
            var e = S(t);
            return e._chain = !0, e
        };
        var $ = function(t, e) {
            return t._chain ? S(e).chain() : e
        };
        S.mixin = function(t) {
            S.each(S.functions(t), function(e) {
                var n = S[e] = t[e];
                S.prototype[e] = function() {
                    var t = [this._wrapped];
                    return p.apply(t, arguments), $(this, n.apply(S, t))
                }
            })
        }, S.mixin(S), S.each(["pop", "push", "reverse", "shift", "sort", "splice", "unshift"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                var n = this._wrapped;
                return e.apply(n, arguments), "shift" !== t && "splice" !== t || 0 !== n.length || delete n[0], $(this, n)
            }
        }), S.each(["concat", "join", "slice"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                return $(this, e.apply(this._wrapped, arguments))
            }
        }), S.prototype.value = function() {
            return this._wrapped
        }, S.prototype.valueOf = S.prototype.toJSON = S.prototype.value, S.prototype.toString = function() {
            return "" + this._wrapped
        }, r = [], i = function() {
            return S
        }.apply(e, r), !(void 0 !== i && (t.exports = i))
    }).call(void 0)
}, function(t, e) {
    "use strict";
    t.exports = {
        render_x_axis_label: function() {
            var t = "Model Position";
            this.display_ali_map && (t = "Alignment Column"), this.called_on.find(".logo_xaxis").remove(), this.called_on.prepend('<div class="logo_xaxis" class="centered" style="margin-left:40px"><p class="xaxis_text" style="width:10em;margin:1em auto">' + t + "</p></div>")
        },
        render_y_axis_label: function() {
            this.dom_element.parent().before('<canvas class="logo_yaxis" height="' + this.options.height + '" width="55"></canvas>');
            var t = this.called_on.find(".logo_yaxis"),
                e = (Math.abs(this.data.max_height), isNaN(this.data.min_height_obs) ? 0 : parseInt(this.data.min_height_obs, 10), null),
                n = "Information Content (bits)";
            e = t[0].getContext("2d"), e.beginPath(), e.moveTo(55, 1), e.lineTo(40, 1), e.moveTo(55, this.info_content_height), e.lineTo(40, this.info_content_height), e.moveTo(55, this.info_content_height / 2), e.lineTo(40, this.info_content_height / 2), e.lineWidth = 1, e.strokeStyle = "#666666", e.stroke(), e.fillStyle = "#666666", e.textAlign = "right", e.font = "bold 10px Arial", e.textBaseline = "top", e.fillText(parseFloat(this.data.max_height).toFixed(1), 38, 0), e.textBaseline = "middle", e.fillText(parseFloat(this.data.max_height / 2).toFixed(1), 38, this.info_content_height / 2), e.fillText("0", 38, this.info_content_height), "score" === this.data.height_calc && (n = "Score (bits)"), e.save(), e.translate(5, this.height / 2 - 20), e.rotate(-Math.PI / 2), e.textAlign = "center", e.font = "normal 12px Arial", e.fillText(n, 1, 0), e.restore(), e.fillText("occupancy", 55, this.info_content_height + 7), this.show_inserts && (e.fillText("ins. prob.", 50, 280), e.fillText("ins. len.", 46, 296))
        }
    }
}, function(t, e) {
    "use strict";
    var n = null;
    t.exports = function() {
        if (!n) {
            var t = document.createElement("canvas");
            n = !(!t.getContext || !t.getContext("2d"))
        }
        return n
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#FF9966",
        C: "#009999",
        D: "#FF0000",
        E: "#CC0033",
        F: "#00FF00",
        G: "#f2f20c",
        H: "#660033",
        I: "#CC9933",
        K: "#663300",
        L: "#FF9933",
        M: "#CC99CC",
        N: "#336666",
        P: "#0099FF",
        Q: "#6666CC",
        R: "#990000",
        S: "#0000FF",
        T: "#00FFFF",
        V: "#FFCC33",
        W: "#66CC66",
        Y: "#006600"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#cbf751",
        C: "#5ec0cc",
        G: "#ffdf59",
        T: "#b51f16",
        U: "#b51f16"
    }
}, function(t, e, n) {
    "use strict";
    var r = n(5);
    t.exports = function(t, e, n) {
        t.find(".logo_settings_switch, .logo_settings .close").on("click", function(t) {
            t.preventDefault(), r(".logo_settings").toggle()
        }), t.find(".logo_reset").on("click", function(t) {
            t.preventDefault(), e.changeZoom({
                target: e.default_zoom
            })
        }), t.find(".logo_change").on("click", function(t) {
            t.preventDefault()
        }), t.find(".logo_zoomin").on("click", function(t) {
            t.preventDefault(), e.changeZoom({
                distance: .1,
                direction: "+"
            })
        }), t.find(".logo_zoomout").on("click", function(t) {
            t.preventDefault(), e.changeZoom({
                distance: .1,
                direction: "-"
            })
        }), t.find(".logo_scale").on("change", function(t) {
            e.toggleScale(this.value)
        }), t.find(".logo_color").on("change", function(t) {
            e.toggleColorscheme(this.value)
        }), t.find(".logo_ali_map").on("change", function(t) {
            e.toggleAliMap(this.value)
        }), t.find(".logo_position").on("change", function() {
            this.value.match(/^\d+$/m) && e.scrollToColumn(this.value, 1)
        }), n.on("dblclick", function(n) {
            console.log("dblclick", e), offset = e.logo_graphic.offset(), x = parseInt(n.pageX - offset.left, 10), window_position = n.pageX - t.parent().offset().left, col = e.columnFromCoordinates(x), console.log("col", col), current = e.zoom, current < 1 ? e.changeZoom({
                target: 1,
                offset: window_position,
                column: col
            }) : e.changeZoom({
                target: .3,
                offset: window_position,
                column: col
            })
        }), r(document).on(t.attr("id") + ".scrolledTo", function(t, n, r, i) {
            e.render({
                target: n
            })
        }), r(document).on("keydown", function(t) {
            t.ctrlKey || (61 !== t.which && 107 !== t.which || (zoom += .1, e.changeZoom({
                distance: .1,
                direction: "+"
            })), 109 !== t.which && 0 !== t.which || (zoom -= .1, e.changeZoom({
                distance: .1,
                direction: "-"
            })))
        })
    }
}, function(t, e, n) {
    "use strict";
    _ = n(75);
    var r = n(77),
        i = n(84),
        o = n(83),
        s = n(2),
        u = n(76),
        a = n(80),
        l = n(82),
        c = n(5);
    t.exports = s.extend({
        options: {
            xaxis: !0,
            yaxis: !0,
            height: 300,
            column_width: 34,
            debug: !0,
            scale_height_enabled: !0,
            scaled_max: !0,
            zoom_buttons: !0,
            colorscheme: "default",
            data: void 0,
            start: 1,
            end: void 0,
            zoom: .4,
            colors: void 0,
            divider: !1,
            show_probs: !1,
            divider_step: 5,
            show_divider: !1,
            border: !1,
            settings: !1,
            scroller: !0,
            positionMarker: !0
        },
        loadDefault: function(t) {
            this.data = t.data, this.display_ali_map = 0, this.alphabet = t.data.alphabet || "dna", this.start = t.start, this.zoom = parseFloat(t.zoom) || .4, this.default_zoom = this.zoom, this.column_width = t.column_width, this.height = t.height, this.canvas_width = 5e3, this.scale_height_enabled = t.scale_height_enabled, this.scrollme = null, this.previous_target = 0, this.rendered = [], this.previous_zoom = 0, void 0 == this.data.max_height && (this.data.max_height = this.calcMaxHeight(this.data.heightArr)), this.data.insert_probs && this.data.delete_probs || (this.options.show_probs = !1), t.scaled_max ? this.data.max_height = t.data.max_height_obs || this.data.max_height || 2 : this.data.max_height = t.data.max_height_theory || this.data.max_height || 2, t.colors ? this.changeColors(t.colors) : "aa" === this.alphabet ? (this.aa_colors = n(78), this.changeColors(this.aa_colors)) : (this.dna_colors = n(79), this.changeColors(this.dna_colors))
        },
        initialize: function(t) {
            if (!r()) return void(this.el.textContent = "Your browser doesn't support canvas.");
            void 0 == t.data && (this.el.textContent = "No data added."), _.extend(this.options, t);
            var e = this.options;
            if (this.loadDefault(e), this.options.show_probs ? this.data.processing && /^observed|weighted/.test(this.data.processing) ? (this.show_inserts = 0, this.info_content_height = this.height - 14) : (this.show_inserts = 1, this.info_content_height = this.height - 44) : this.info_content_height = this.height, this.$el = c(this.el), this.initDivs(), this.options.settings) {
                var n = l(this, e);
                this.$el.append(n)
            }
            a(this.$el, this, this.logo_graphic)
        },
        initDivs: function() {
            var t = f("div");
            t.className = "logo_graphic", this.logo_graphic = c(t);
            var e = f("div");
            if (e.className = "logo_container", e.style.height = this.height, this.container = c(e), this.container.append(t), this.$el.append(e), this.options.divider) {
                var n = f("div");
                n.className = "logo_divider", this.$el.append(n)
            }
            this.dom_element = c(t), this.called_on = this.$el, this.options.xaxis && u.render_x_axis_label.call(this), this.options.yaxis ? u.render_y_axis_label.call(this) : this.container[0].style.marginLeft = "0px"
        },
        render: function() {
            return i.call(this), this
        },
        changeColors: function(t) {
            this.colors = t, void 0 != t && void 0 != t.type && (this.colorscheme = "dynamic"), this.buildAlphabet()
        },
        buildAlphabet: function() {
            this.letters = {};
            var t = this.colors;
            if ("dynamic" == this.colorscheme) {
                var e = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
                t = {}, e.forEach(function(e) {
                    t[e] = ""
                })
            }
            for (var n in t)
                if (t.hasOwnProperty(n)) {
                    var r = {
                        color: t[n]
                    };
                    this.letters[n] = new o(n, r)
                }
        },
        toggleColorscheme: function(t) {
            var e = this.currentColumn();
            t ? "default" === t ? this.colorscheme = "default" : this.colorscheme = "consensus" : "default" === this.colorscheme ? this.colorscheme = "consensus" : this.colorscheme = "default", this.rendered = [], this.scrollme.reflow(), this.scrollToColumn(e + 1), this.scrollToColumn(e)
        },
        toggleScale: function(t) {
            var e = this.currentColumn();
            t ? "obs" === t ? this.data.max_height = this.data.max_height_obs : this.data.max_height = this.data.max_height_theory : this.data.max_height === this.data.max_height_obs ? this.data.max_height = this.data.max_height_theory : this.data.max_height = this.data.max_height_obs, this.rendered = [], this.logoYAxis && this.logoYAxis.remove(), u.render_y_axis_label.call(this), this.scrollme.reflow(), this.scrollToColumn(e + 1), this.scrollToColumn(e)
        },
        toggleAliMap: function(t) {
            var e = this.currentColumn();
            t ? "model" === t ? this.display_ali_map = 0 : this.display_ali_map = 1 : 1 === this.display_ali_map ? this.display_ali_map = 0 : this.display_ali_map = 1, u.render_x_axis_label(this), this.rendered = [], this.scrollme.reflow(), this.scrollToColumn(e + 1), this.scrollToColumn(e)
        },
        currentColumn: function() {
            var t = this.scrollme.scroller.getValues().left,
                e = this.column_width * this.zoom,
                n = t / e,
                r = this.container.width() / e / 2;
            return Math.ceil(n + r)
        },
        changeZoom: function(t) {
            var e = .3,
                n = null;
            if (t.target ? e = t.target : t.distance && (e = (parseFloat(this.zoom) - parseFloat(t.distance)).toFixed(1), "+" === t.direction && (e = (parseFloat(this.zoom) + parseFloat(t.distance)).toFixed(1))), e > 1 ? e = 1 : e < .1 && (e = .1), n = this.logo_graphic.width() * e / this.zoom, n > this.container.width())
                if (t.column) {
                    this.zoom = e, this.render({
                        zoom: this.zoom
                    }), this.scrollme.reflow();
                    var r = this.coordinatesFromColumn(t.column);
                    this.scrollme.scroller.scrollTo(r - t.offset)
                } else {
                    var i = this.currentColumn();
                    this.zoom = e, this.render({
                        zoom: this.zoom
                    }), this.scrollme.reflow(), this.scrollToColumn(i)
                } return this.zoom
        },
        columnFromCoordinates: function(t) {
            return Math.ceil(t / (this.column_width * this.zoom))
        },
        coordinatesFromColumn: function(t) {
            return (t - 1) * (this.column_width * this.zoom) + this.column_width * this.zoom / 2
        },
        scrollToColumn: function(t, e) {
            var n = this.logo_container.width() / 2,
                r = this.coordinatesFromColumn(t);
            this.scrollme.scroller.scrollTo(r - n, 0, e)
        },
        calcMaxHeight: function(t) {
            return t.reduce(function(t, e) {
                var n = 0;
                for (var r in e) n += e[r];
                return n > t ? n : t
            }, 0)
        }
    });
    var f = function(t) {
        return document.createElement(t)
    }
}, function(t, e, n) {
    "use strict";
    var r = n(5);
    t.exports = function(t, e) {
        var n = r('<form class="logo_form"><fieldset><label for="position">Column number</label><input type="text" name="position" class="logo_position"></input><button class="button logo_change">Go</button></fieldset></form>'),
            i = r('<div class="logo_settings"></div>');
        if (i.append('<span class="close">x</span>'), t.scale_height_enabled && t.data.max_height_obs < t.data.max_height_theory) {
            var o = "",
                s = "",
                u = "",
                a = "";
            t.data.max_height_obs === t.data.max_height ? o = "checked" : s = "checked"
        }
        var l = '<fieldset><legend>Scale</legend><label><input type="radio" name="scale" class="logo_scale" value="obs" ' + o + "/>Maximum Observed " + a + '</label></br><label><input type="radio" name="scale" class="logo_scale" value="theory" ' + s + "/>Maximum Theoretical " + u + "</label></fieldset>";
        if (i.append(l), "score" !== t.data.height_calc && "aa" === t.data.alphabet && t.data.probs_arr) {
            var c = null,
                f = null,
                h = "",
                d = "";
            "default" === t.colorscheme ? c = "checked" : f = "checked", e.help && (h = '<a class="help" href="/help#colors_default" title="Each letter receives its own color."><span aria-hidden="true" data-icon="?"></span><span class="reader-text">help</span></a>', d = '<a class="help" href="/help#colors_consensus" title="Letters are colored as in Clustalx and Jalview, with colors depending on composition of the column."><span aria-hidden="true" data-icon="?"></span><span class="reader-text">help</span></a>');
            var p = '<fieldset><legend>Color Scheme</legend><label><input type="radio" name="color" class="logo_color" value="default" ' + c + "/>Default " + h + '</label></br><label><input type="radio" name="color" class="logo_color" value="consensus" ' + f + "/>Consensus Colors " + d + "</label></fieldset>";
            i.append(p)
        }
        if (t.data.ali_map) {
            var g = null,
                v = null,
                m = "",
                y = "";
            0 === t.display_ali_map ? g = "checked" : v = "checked", e.help && (m = '<a class="help" href="/help#coords_model" title="The coordinates along the top of the plot show the model position."><span aria-hidden="true" data-icon="?"></span><span class="reader-text">help</span></a>', y = '<a class="help" href="/help#coords_ali" title="The coordinates along the top of the plot show the column in the alignment associated with the model"><span aria-hidden="true" data-icon="?"></span><span class="reader-text">help</span></a>');
            var _ = '<fieldset><legend>Coordinates</legend><label><input type="radio" name="coords" class="logo_ali_map" value="model" ' + g + "/>Model " + m + '</label></br><label><input type="radio" name="coords" class="logo_ali_map" value="alignment" ' + v + "/>Alignment " + y + "</label></fieldset>";
            i.append(_)
        }
        var b = r('<div class="logo_controls"></div>');
        return t.zoom_enabled && b.append('<button class="logo_zoomout button">-</button><button class="logo_zoomin button">+</button>'), i.children().length > 0 && (b.append('<button class="logo_settings_switch button">Settings</button>'), b.append(i)), n.append(b), n
    }
}, function(t, e) {
    "use strict";
    t.exports = function(t, e) {
        e = e || {}, this.value = t, this.width = parseInt(e.width, 10) || 100, "W" === this.value && (this.width += 30 * this.width / 100), this.height = parseInt(e.height, 10) || 100, this.color = e.color || "#000000", this.fontSize = e.fontSize || 138, this.scaled = function() {}, this.draw = function(t, e, n, r, i, o) {
            var s = e / this.height,
                u = n / this.width,
                a = t.font;
            t.transform(u, 0, 0, s, r, i), t.fillStyle = o || this.color, t.textAlign = "center", t.font = "bold " + this.fontSize + "px Arial", t.fillText(this.value, 0, 0), t.setTransform(1, 0, 0, 1, 0, 0), t.fillStyle = "#000000", t.font = a
        }
    }
}, function(t, e, n) {
    "use strict";

    function r(t, e, n, r, i) {
        var o = s(t).find("#canv_" + r);
        return o.length || (s(t).append('<canvas class="canvas_logo" id="canv_' + r + '"  height="' + e + '" width="' + n + '" style="left:' + i * r + 'px"></canvas>'),
            o = s(t).find("#canv_" + r)), s(o).attr("width", n).attr("height", e), o[0]
    }
    var i = n(86),
        o = n(85),
        s = n(5);
    t.exports = function(t) {
        if (this.data) {
            t = t || {};
            var e = t.zoom || this.zoom,
                n = t.target || 1,
                s = (t.scaled || null, this.dom_element.parent().attr("width")),
                u = 1,
                a = null,
                l = null,
                c = 0;
            if (this.previous_target = n, t.start && (this.start = t.start), t.end && (this.end = t.end), e <= .1 ? e = .1 : e >= 1 && (e = 1), this.zoom = e, a = this.end || this.data.heightArr.length, l = this.start || 1, a = a > this.data.heightArr.length ? this.data.heightArr.length : a, a = a < l ? l : a, l = l > a ? a : l, l = l > 1 ? l : 1, this.y = this.height - 20, this.max_width = this.column_width * (a - l + 1), s > this.max_width && (e = 1, this.zoom_enabled = !1), this.zoom = e, this.zoomed_column = this.column_width * e, this.total_width = this.zoomed_column * (a - l + 1), e < 1)
                for (; this.total_width < s && (this.zoom += .1, this.zoomed_column = this.column_width * this.zoom, this.total_width = this.zoomed_column * (a - l + 1), this.zoom_enabled = !1, !(e >= 1)););
            n > this.total_width && (n = this.total_width), this.dom_element.attr({
                width: this.total_width + "px"
            }).css({
                width: this.total_width + "px"
            }), this.canvas_width = this.total_width;
            var f = Math.ceil(this.total_width / this.canvas_width);
            for (this.columns_per_canvas = Math.ceil(this.canvas_width / this.zoomed_column), this.previous_zoom !== this.zoom && (this.dom_element.find("canvas").remove(), this.previous_zoom = this.zoom, this.rendered = []), this.canvases = [], this.contexts = [], c = 0; c < f; c++) {
                var h = this.columns_per_canvas * c + l,
                    d = h + this.columns_per_canvas - 1;
                d > a && (d = a);
                var p = (d - h + 1) * this.zoomed_column;
                p > u && (u = p);
                var g = u * c,
                    v = g + p;
                if (n < v + v / 2 && n > g - g / 2)
                    if (this.canvases[c] = r(this.dom_element, this.height, p, c, u), this.contexts[c] = this.canvases[c].getContext("2d"), this.contexts[c].setTransform(1, 0, 0, 1, 0, 0), this.contexts[c].clearRect(0, 0, p, this.height), this.contexts[c].fillStyle = "#ffffff", this.contexts[c].fillRect(0, 0, v, this.height), this.zoomed_column > 12) {
                        var m = parseInt(10 * e, 10);
                        m = m > 10 ? 10 : m, this.debug && o.call(this, h, d, c, 1), i.call(this, h, d, c, m)
                    } else o.call(this, h, d, c)
            }!this.scrollme && this.options.scroller && (this.scrollme = new EasyScroller(this.dom_element[0], {
                scrollingX: 1,
                scrollingY: 0,
                eventTarget: this.called_on
            })), 1 !== n && this.scrollme.reflow()
        }
    }
}, function(t, e, n) {
    "use strict";

    function r(t, e, n, r, i, o, s, u) {
        var a = "#ffffff";
        u ? (i > .1 ? a = "#d7301f" : i > .05 ? a = "#fc8d59" : i > .03 && (a = "#fdcc8a"), t.fillStyle = a, t.fillRect(e, n + 15, r, 10), a = "#ffffff", o > 9 ? a = "#d7301f" : o > 7 ? a = "#fc8d59" : o > 4 && (a = "#fdcc8a"), t.fillStyle = a, t.fillRect(e, n + 30, r, 10)) : n += 30, a = "#ffffff", s < .75 ? a = "#2171b5" : s < .85 ? a = "#6baed6" : s < .95 && (a = "#bdd7e7"), t.fillStyle = a, t.fillRect(e, n, r, 10)
    }
    var i = n(44),
        o = n(46),
        s = n(45);
    t.exports = function(t, e, n, u) {
        var a = 0,
            l = t,
            c = null,
            f = 0,
            h = Math.abs(this.data.max_height),
            d = Math.abs(this.data.min_height_obs),
            p = h + d,
            g = Math.round(100 * Math.abs(this.data.max_height) / p),
            v = Math.round(this.info_content_height * g / 100),
            m = (this.info_content_height - v, 10);
        for (f = t; f <= e; f++) {
            if (this.data.mmline && 1 === this.data.mmline[f - 1]) this.contexts[n].fillStyle = "#cccccc", this.contexts[n].fillRect(a, 10, this.zoomed_column, this.height - 40);
            else {
                var y = this.data.heightArr[f - 1],
                    _ = 0,
                    b = (y.length, 0);
                for (var b in y) {
                    var x = [b, y[b]];
                    if (x[1] > .01) {
                        var w = parseFloat(x[1]) / this.data.max_height,
                            S = a,
                            k = (this.info_content_height - 2) * w,
                            j = this.info_content_height - 2 - _ - k,
                            O = null;
                        O = "dynamic" === this.colorscheme ? this.colors.getColor(x[0], {
                            pos: f - 1
                        }) : "consensus" === this.colorscheme ? this.cmap[f - 1][x[0]] || "#7a7a7a" : this.colors[x[0]], u ? (this.contexts[n].strokeStyle = O, this.contexts[n].strokeRect(S, j, this.zoomed_column, k)) : (this.contexts[n].fillStyle = O, this.contexts[n].fillRect(S, j, this.zoomed_column, k)), _ += k
                    }
                }
            }
            this.zoom < .2 ? m = 20 : this.zoom < .3 && (m = 10), this.options.positionMarker && f % m === 0 && (this.options.show_probs && o(this.contexts[n], a + this.zoomed_column, this.height - 30, parseFloat(this.height), "#dddddd"), o(this.contexts[n], a + this.zoomed_column, 0, 5), c = this.display_ali_map ? this.data.ali_map[f - 1] : l, s(this.contexts[n], a - 2, 10, this.zoomed_column, c, 10, !0)), this.options.show_probs && r(this.contexts[n], a, this.height - 42, this.zoomed_column, this.data.insert_probs[f - 1], this.data.insert_lengths[f - 1], this.data.delete_probs[f - 1], this.show_inserts), this.options.show_probs && (this.show_inserts ? i(this.contexts[n], this.height - 45, this.total_width) : i(this.contexts[n], this.height - 15, this.total_width)), this.options.border && i(this.contexts[n], 0, this.total_width), a += this.zoomed_column, l++
        }
    }
}, function(t, e, n) {
    "use strict";

    function r(t, e, n, r, o, s, u) {
        var a = n - 4,
            l = "#ffffff",
            c = "#555555";
        u && (a = n - 35), o < .75 ? (l = "#2171b5", c = "#ffffff") : o < .85 ? l = "#6baed6" : o < .95 && (l = "#bdd7e7"), i(t, e, a, o, s, r, l, c)
    }

    function i(t, e, n, r, i, o, s, u) {
        t.font = i + "px Arial", t.fillStyle = s, t.fillRect(e, n - 10, o, 14), t.textAlign = "center", t.fillStyle = u, t.fillText(r, e + o / 2, n)
    }

    function o(t, e) {
        var n = e.ralign ? e.x + t.zoomed_column : e.x,
            r = e.ralign ? e.x + 2 : e.x;
        l(t.contexts[e.context_num], n, t.height - 30, -30 - t.height, "#dddddd"), l(t.contexts[e.context_num], n, 0, 5), c(t.contexts[e.context_num], r, 10, t.zoomed_column, e.column_num, e.fontsize, e.ralign)
    }

    function s(t, e, n, r, o, s) {
        var u = n - 20,
            a = "#ffffff",
            c = "#555555";
        o > .1 ? (a = "#d7301f", c = "#ffffff") : o > .05 ? a = "#fc8d59" : o > .03 && (a = "#fdcc8a"), i(t, e, u, o, s, r, a, c), o > .03 && l(t, e + r, n - 30, -30 - n, a)
    }

    function u(t, e, n, r, o, s) {
        var u = "#ffffff",
            a = "#555555";
        o > 9 ? (u = "#d7301f", a = "#ffffff") : o > 7 ? u = "#fc8d59" : o > 4 && (u = "#fdcc8a"), i(t, e, n, o, s, r, u, a)
    }
    var a = n(44),
        l = n(46),
        c = n(45);
    t.exports = function(t, e, n, i) {
        var c = 0,
            f = t,
            h = null,
            d = 0,
            p = Math.abs(this.data.max_height),
            g = isNaN(this.data.min_height_obs) ? 0 : parseInt(this.data.min_height_obs, 10),
            v = p + Math.abs(g),
            m = Math.round(100 * Math.abs(this.data.max_height) / v),
            y = Math.round(this.info_content_height * m / 100),
            _ = this.info_content_height - y;
        for (y / this.info_content_height, _ / this.info_content_height, e + 3 <= this.end && (e += 3), d = t; d <= e; d++) {
            if (this.data.mmline && 1 === this.data.mmline[d - 1]) this.contexts[n].fillStyle = "#cccccc", this.contexts[n].fillRect(c, 10, this.zoomed_column, this.height - 40);
            else {
                var b = this.data.heightArr[d - 1],
                    x = [];
                if (b) {
                    var w = 0,
                        S = (b.length, 0),
                        k = null;
                    for (var S in b) {
                        var j = b[S],
                            O = [S, j],
                            E = c + this.zoomed_column / 2,
                            M = null;
                        if (O[1] > .01) {
                            M = parseFloat(O[1]) / this.data.max_height;
                            var z = this.info_content_height - 2 - w,
                                A = (this.info_content_height - 2) * M;
                            x[S] = [A, this.zoomed_column, E, z], w += A
                        }
                    }
                    for (var S in b) x[S] && this.letters[S] && (k = "dynamic" === this.colorscheme ? this.colors.getColor(S, {
                        pos: d - 1
                    }) : "consensus" === this.colorscheme ? this.cmap[d - 1][S] || "#7a7a7a" : null, this.letters[S].draw(this.contexts[n], x[S][0], x[S][1], x[S][2], x[S][3], k))
                }
            }
            h = this.display_ali_map ? this.data.ali_map[d - 1] : f, this.options.show_divider && (this.zoom < .7 ? d % this.options.divider_step === 0 && o(this, {
                context_num: n,
                x: c,
                fontsize: 10,
                column_num: h,
                ralign: !0
            }) : o(this, {
                context_num: n,
                x: c,
                fontsize: i,
                column_num: h
            })), this.options.show_probs && (r(this.contexts[n], c, this.height, this.zoomed_column, this.data.delete_probs[d - 1], i, this.show_inserts), l(this.contexts[n], c, this.height - 15, 5), this.show_inserts && (s(this.contexts[n], c, this.height, this.zoomed_column, this.data.insert_probs[d - 1], i), u(this.contexts[n], c, this.height - 5, this.zoomed_column, this.data.insert_lengths[d - 1], i), l(this.contexts[n], c, this.height - 45, 5), l(this.contexts[n], c, this.height - 30, 5))), c += this.zoomed_column, f++
        }
        this.options.show_probs && (this.show_inserts && (a(this.contexts[n], this.height - 30, this.total_width), a(this.contexts[n], this.height - 45, this.total_width)), a(this.contexts[n], this.height - 15, this.total_width)), this.options.border && a(this.contexts[n], 0, this.total_width)
    }
}, function(t, e) {
    "use strict";
    var n = window.HTMLCanvasElement && window.HTMLCanvasElement.prototype,
        r = window.Blob && function() {
            try {
                return Boolean(new Blob)
            } catch (t) {
                return !1
            }
        }(),
        i = r && window.Uint8Array && function() {
            try {
                return 100 === new Blob([new Uint8Array(100)]).size
            } catch (t) {
                return !1
            }
        }(),
        o = window.BlobBuilder || window.WebKitBlobBuilder || window.MozBlobBuilder || window.MSBlobBuilder,
        s = (r || o) && window.atob && window.ArrayBuffer && window.Uint8Array && function(t) {
            var e, n, s, u, a, l;
            for (e = t.split(",")[0].indexOf("base64") >= 0 ? atob(t.split(",")[1]) : decodeURIComponent(t.split(",")[1]), n = new ArrayBuffer(e.length), s = new Uint8Array(n), u = 0; u < e.length; u += 1) s[u] = e.charCodeAt(u);
            return a = t.split(",")[0].split(":")[1].split(";")[0], r ? new Blob([i ? s : n], {
                type: a
            }) : (l = new o, l.append(n), l.getBlob(a))
        };
    window.HTMLCanvasElement && !n.toBlob && (n.mozGetAsFile ? n.toBlob = function(t, e, r) {
        t(r && n.toDataURL && s ? s(this.toDataURL(e, r)) : this.mozGetAsFile("blob", e))
    } : n.toDataURL && s && (n.toBlob = function(t, e, n) {
        t(s(this.toDataURL(e, n)))
    })), t.exports = s
}, function(t, e, n) {
    (function(t) {
        "use strict"; /*! @source http://purl.eligrey.com/github/FileSaver.js/blob/master/FileSaver.js */
        var e = e || "undefined" != typeof navigator && navigator.msSaveOrOpenBlob && navigator.msSaveOrOpenBlob.bind(navigator) || function(t) {
                if ("undefined" == typeof navigator || !/MSIE [1-9]\./.test(navigator.userAgent)) {
                    var e = t.document,
                        n = function() {
                            return t.URL || t.webkitURL || t
                        },
                        r = e.createElementNS("http://www.w3.org/1999/xhtml", "a"),
                        i = !t.externalHost && "download" in r,
                        o = function(n) {
                            var r = e.createEvent("MouseEvents");
                            r.initMouseEvent("click", !0, !1, t, 0, 0, 0, 0, 0, !1, !1, !1, !1, 0, null), n.dispatchEvent(r)
                        },
                        s = t.webkitRequestFileSystem,
                        u = t.requestFileSystem || s || t.mozRequestFileSystem,
                        a = function(e) {
                            (t.setImmediate || t.setTimeout)(function() {
                                throw e
                            }, 0)
                        },
                        l = "application/octet-stream",
                        c = 0,
                        f = [],
                        h = function() {
                            for (var t = f.length; t--;) {
                                var e = f[t];
                                "string" == typeof e ? n().revokeObjectURL(e) : e.remove()
                            }
                            f.length = 0
                        },
                        d = function(t, e, n) {
                            e = [].concat(e);
                            for (var r = e.length; r--;) {
                                var i = t["on" + e[r]];
                                if ("function" == typeof i) try {
                                    i.call(t, n || t)
                                } catch (o) {
                                    a(o)
                                }
                            }
                        },
                        p = function(e, a) {
                            var h, p, g, v = this,
                                m = e.type,
                                y = !1,
                                _ = function() {
                                    var t = n().createObjectURL(e);
                                    return f.push(t), t
                                },
                                b = function() {
                                    d(v, "writestart progress write writeend".split(" "))
                                },
                                x = function() {
                                    !y && h || (h = _(e)), p ? p.location.href = h : window.open(h, "_blank"), v.readyState = v.DONE, b()
                                },
                                w = function(t) {
                                    return function() {
                                        if (v.readyState !== v.DONE) return t.apply(this, arguments)
                                    }
                                },
                                S = {
                                    create: !0,
                                    exclusive: !1
                                };
                            return v.readyState = v.INIT, a || (a = "download"), i ? (h = _(e), r.href = h, r.download = a, o(r), v.readyState = v.DONE, void b()) : (t.chrome && m && m !== l && (g = e.slice || e.webkitSlice, e = g.call(e, 0, e.size, l), y = !0), s && "download" !== a && (a += ".download"), (m === l || s) && (p = t), u ? (c += e.size, void u(t.TEMPORARY, c, w(function(t) {
                                t.root.getDirectory("saved", S, w(function(t) {
                                    var n = function() {
                                        t.getFile(a, S, w(function(t) {
                                            t.createWriter(w(function(n) {
                                                n.onwriteend = function(e) {
                                                    p.location.href = t.toURL(), f.push(t), v.readyState = v.DONE, d(v, "writeend", e)
                                                }, n.onerror = function() {
                                                    var t = n.error;
                                                    t.code !== t.ABORT_ERR && x()
                                                }, "writestart progress write abort".split(" ").forEach(function(t) {
                                                    n["on" + t] = v["on" + t]
                                                }), n.write(e), v.abort = function() {
                                                    n.abort(), v.readyState = v.DONE
                                                }, v.readyState = v.WRITING
                                            }), x)
                                        }), x)
                                    };
                                    t.getFile(a, {
                                        create: !1
                                    }, w(function(t) {
                                        t.remove(), n()
                                    }), w(function(t) {
                                        t.code === t.NOT_FOUND_ERR ? n() : x()
                                    }))
                                }), x)
                            }), x)) : void x())
                        },
                        g = p.prototype,
                        v = function(t, e) {
                            return new p(t, e)
                        };
                    return g.abort = function() {
                        var t = this;
                        t.readyState = t.DONE, d(t, "abort")
                    }, g.readyState = g.INIT = 0, g.WRITING = 1, g.DONE = 2, g.error = g.onwritestart = g.onprogress = g.onwrite = g.onabort = g.onerror = g.onwriteend = null, t.addEventListener("unload", h, !1), v.unload = function() {
                        h(), t.removeEventListener("unload", h, !1)
                    }, v
                }
            }("undefined" != typeof self && self || "undefined" != typeof window && window || (void 0).content),
            n = window.define;
        "undefined" == typeof n && "undefined" != typeof window.almond && "define" in window.almond && (n = window.almond.define), "undefined" != typeof t && null !== t ? t.exports = e : "undefined" != typeof n && null !== n && null != n.amd && n("saveAs", [], function() {
            return e
        })
    }).call(e, n(20)(t))
}, function(t, e) {
    "use strict";
    t.exports = function() {
        var t = [];
        return t.toString = function() {
            for (var t = [], e = 0; e < this.length; e++) {
                var n = this[e];
                n[2] ? t.push("@media " + n[2] + "{" + n[1] + "}") : t.push(n[1])
            }
            return t.join("")
        }, t.i = function(e, n) {
            "string" == typeof e && (e = [
                [null, e, ""]
            ]);
            for (var r = {}, i = 0; i < this.length; i++) {
                var o = this[i][0];
                "number" == typeof o && (r[o] = !0)
            }
            for (i = 0; i < e.length; i++) {
                var s = e[i];
                "number" == typeof s[0] && r[s[0]] || (n && !s[2] ? s[2] = n : n && (s[2] = "(" + s[2] + ") and (" + n + ")"), t.push(s))
            }
        }, t
    }
}, function(t, e, n) {
    "use strict";
    ! function(n) {
        function r(t, e) {
            if (!(this instanceof r)) return new r(t, e);
            this.domain = [], this.range = [], Array.isArray(t) && (this.domain = t), Array.isArray(e) && (this.range = e);
            var n = function(t) {
                if ("number" != typeof t) return null;
                var e = this.domain[0],
                    n = this.domain[1],
                    r = this.range[0],
                    i = this.range[1];
                return "number" !== r && "number" != typeof i && (r = e, i = n), r + (i - r) / (n - e) * (t - e)
            }.bind(this);
            return n.domain = function(t) {
                return Array.isArray(t) && (this.domain = t), n
            }.bind(this), n.range = function(t) {
                return Array.isArray(t) && (this.range = t), n
            }.bind(this), n
        }
        "undefined" != typeof t && t.exports && (e = t.exports = r), e.LinearScale = r
    }(void 0)
}, function(t, e, n) {
    "use strict";
    t.exports = n(92)
}, function(t, e, n) {
    "use strict";
    var r, i, o;
    i = n(5), o = n(2), t.exports = r = o.extend({
        initialize: function(t) {
            this._nodes = [], this.name = t.name || "", this.el.className += "smenubar"
        },
        render: function() {
            for (var t = this.el.firstChild; t;) this.el.removeChild(t), t = this.el.firstChild;
            this.el.appendChild(this.buildDOM())
        },
        setName: function(t) {
            this.name = t
        },
        addNode: function(t, e, n) {
            var r;
            null != n && (r = n.style), null == this._nodes && (this._nodes = []), this._nodes.push({
                label: t,
                callback: e,
                style: r
            })
        },
        getNode: function(t) {
            var e = void 0;
            return this._nodes.forEach(function(n) {
                n.label === t && (e = n)
            }), e
        },
        modifyNode: function(t, e, n) {
            var r = this.getNode(t);
            r.callback = e || r.callback, n = n || {}, r.style = n.style || r.style
        },
        renameNode: function(t, e) {
            var n = this.getNode(t);
            n.label = e || n.label
        },
        removeNode: function(t) {
            var e = this.getNode(t);
            this._nodes.splice(this._nodes.indexOf(e), 1)
        },
        removeAllNodes: function() {
            this._nodes = []
        },
        buildDOM: function() {
            var t = document.createElement("span");
            return t.appendChild(this._buildM({
                nodes: this._nodes,
                name: this.name
            })), t
        },
        _buildM: function(t) {
            var e, n, r, o, s, u, a, l = t.nodes,
                c = t.name,
                f = document.createElement("div");
            f.className = "smenu-dropdown smenu-dropdown-tip", f.style.display = "none";
            var h = document.createElement("ul");
            h.className = "smenu-dropdown-menu";
            for (var d = 0, p = l.length; d < p; d++) {
                s = l[d], o = document.createElement("li"), o.textContent = s.label, a = s.style;
                for (r in a) u = a[r], o.style[r] = u;
                o.addEventListener("click", s.callback), this.trigger("new:node", o), h.appendChild(o)
            }
            return this.trigger("new:menu", h), f.appendChild(h), e = document.createElement("a"), e.textContent = c, e.className = "smenubar_alink", this.trigger("new:button", e), i(e).on("click", function(t) {
                return function(n) {
                    return t._showMenu(n, f, e), window.setTimeout(function() {
                        return i(document.body).one("click", function(t) {
                            return f.style.display = "none"
                        })
                    }, 5)
                }
            }(this)), n = document.createDocumentFragment(), n.appendChild(f), n.appendChild(e), n
        },
        _showMenu: function(t, e, n) {
            var r;
            e.style.display = "block", e.style.position = "absolute", r = n.getBoundingClientRect(), e.style.left = r.left + "px", e.style.top = r.top + n.offsetHeight + "px"
        }
    })
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#00a35c",
        R: "#00fc03",
        N: "#00eb14",
        D: "#00eb14",
        C: "#0000ff",
        Q: "#00f10e",
        E: "#00f10e",
        G: "#009d62",
        H: "#00d52a",
        I: "#0054ab",
        L: "#007b84",
        K: "#00ff00",
        M: "#009768",
        F: "#008778",
        P: "#00e01f",
        S: "#00d52a",
        T: "#00db24",
        W: "#00a857",
        Y: "#00e619",
        V: "#005fa0",
        B: "#00eb14",
        X: "#00b649",
        Z: "#00f10e"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#BBBBBB",
        B: "grey",
        C: "yellow",
        D: "red",
        E: "red",
        F: "magenta",
        G: "brown",
        H: "#00FFFF",
        I: "#BBBBBB",
        J: "#fff",
        K: "#00FFFF",
        L: "#BBBBBB",
        M: "#BBBBBB",
        N: "green",
        O: "#fff",
        P: "brown",
        Q: "green",
        R: "#00FFFF",
        S: "green",
        T: "green",
        U: "#fff",
        V: "#BBBBBB",
        W: "magenta",
        X: "grey",
        Y: "magenta",
        Z: "grey",
        Gap: "grey"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "orange",
        B: "#fff",
        C: "green",
        D: "red",
        E: "red",
        F: "blue",
        G: "orange",
        H: "red",
        I: "green",
        J: "#fff",
        K: "red",
        L: "green",
        M: "green",
        N: "#fff",
        O: "#fff",
        P: "orange",
        Q: "#fff",
        R: "red",
        S: "orange",
        T: "orange",
        U: "#fff",
        V: "green",
        W: "blue",
        X: "#fff",
        Y: "blue",
        Z: "#fff",
        Gap: "#fff"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#80a0f0",
        R: "#f01505",
        N: "#00ff00",
        D: "#c048c0",
        C: "#f08080",
        Q: "#00ff00",
        E: "#c048c0",
        G: "#f09048",
        H: "#15a4a4",
        I: "#80a0f0",
        L: "#80a0f0",
        K: "#f01505",
        M: "#80a0f0",
        F: "#80a0f0",
        P: "#ffff00",
        S: "#00ff00",
        T: "#00ff00",
        W: "#80a0f0",
        Y: "#15a4a4",
        V: "#80a0f0",
        B: "#fff",
        X: "#fff",
        Z: "#fff"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#e718e7",
        R: "#6f906f",
        N: "#1be41b",
        D: "#778877",
        C: "#23dc23",
        Q: "#926d92",
        E: "#ff00ff",
        G: "#00ff00",
        H: "#758a75",
        I: "#8a758a",
        L: "#ae51ae",
        K: "#a05fa0",
        M: "#ef10ef",
        F: "#986798",
        P: "#00ff00",
        S: "#36c936",
        T: "#47b847",
        W: "#8a758a",
        Y: "#21de21",
        V: "#857a85",
        B: "#49b649",
        X: "#758a75",
        Z: "#c936c9"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#ad0052",
        B: "#0c00f3",
        C: "#c2003d",
        D: "#0c00f3",
        E: "#0c00f3",
        F: "#cb0034",
        G: "#6a0095",
        H: "#1500ea",
        I: "#ff0000",
        J: "#fff",
        K: "#0000ff",
        L: "#ea0015",
        M: "#b0004f",
        N: "#0c00f3",
        O: "#fff",
        P: "#4600b9",
        Q: "#0c00f3",
        R: "#0000ff",
        S: "#5e00a1",
        T: "#61009e",
        U: "#fff",
        V: "#f60009",
        W: "#5b00a4",
        X: "#680097",
        Y: "#4f00b0",
        Z: "#0c00f3"
    }
}, function(t, e, n) {
    "use strict";

    function r(t) {
        if (null == t || "object" != ("undefined" == typeof t ? "undefined" : i(t))) return t;
        var e = t.constructor();
        for (var n in t) t.hasOwnProperty(n) && (e[n] = t[n]);
        return e
    }
    var i = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        },
        o = n(105),
        s = o.stat,
        u = o.dyn,
        a = n(93),
        l = n(94),
        c = n(95),
        f = n(96),
        h = n(97),
        d = n(98),
        p = n(100),
        g = n(101),
        v = n(102),
        m = n(104),
        y = n(106),
        _ = n(107),
        b = n(108),
        x = n(109),
        w = {
            buried: a,
            buried_index: a,
            cinema: l,
            clustal2: f,
            clustal: c,
            helix: h,
            helix_propensity: h,
            hydro: d,
            lesk: p,
            mae: g,
            nucleotide: v,
            purine: m,
            purine_pyrimidine: m,
            strand: y,
            strand_propensity: y,
            taylor: _,
            turn: b,
            turn_propensity: b,
            zappo: x
        },
        S = n(103),
        k = {
            pid: S
        },
        j = function(t) {
            this.maps = r(w), this.dyn = r(k), this.opt = t
        };
    t.exports = j, j.getScheme = function(t) {
        return w[t]
    }, j.prototype.getScheme = function(t) {
        var e = this.maps[t];
        return void 0 === e && (e = {}, void 0 != this.dyn[t]) ? new u(this.dyn[t], this.opt) : new s(e)
    }, j.prototype.addStaticScheme = function(t, e) {
        this.maps[t] = e
    }, j.prototype.addDynScheme = function(t, e) {
        this.dyn[t] = e
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: " orange",
        B: " #fff",
        C: " green",
        D: " red",
        E: " red",
        F: " green",
        G: " orange",
        H: " magenta",
        I: " green",
        J: " #fff",
        K: " red",
        L: " green",
        M: " green",
        N: " magenta",
        O: " #fff",
        P: " green",
        Q: " magenta",
        R: " red",
        S: " orange",
        T: " orange",
        U: " #fff",
        V: " green",
        W: " green",
        X: " #fff",
        Y: " green",
        Z: " #fff",
        Gap: " #fff"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: " #77dd88",
        B: " #fff",
        C: " #99ee66",
        D: " #55bb33",
        E: " #55bb33",
        F: " #9999ff",
        G: " #77dd88",
        H: " #5555ff",
        I: " #66bbff",
        J: " #fff",
        K: " #ffcc77",
        L: " #66bbff",
        M: " #66bbff",
        N: " #55bb33",
        O: " #fff",
        P: " #eeaaaa",
        Q: " #55bb33",
        R: " #ffcc77",
        S: " #ff4455",
        T: " #ff4455",
        U: " #fff",
        V: " #66bbff",
        W: " #9999ff",
        X: " #fff",
        Y: " #9999ff",
        Z: " #fff",
        Gap: " #fff"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: " #64F73F",
        C: " #FFB340",
        G: " #EB413C",
        T: " #3C88EE",
        U: " #3C88EE"
    }
}, function(t, e) {
    "use strict";
    var n;
    t.exports = n = {}, n.init = function() {
        this.cons = this.opt.conservation()
    }, n.run = function(t, e) {
        var n = this.cons[e.pos];
        return n > .8 ? "#6464ff" : n > .6 ? "#9da5ff" : n > .4 ? "#cccccc" : "#ffffff"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: " #FF83FA",
        C: " #40E0D0",
        G: " #FF83FA",
        R: " #FF83FA",
        T: " #40E0D0",
        U: " #40E0D0",
        Y: " #40E0D0"
    }
}, function(t, e) {
    "use strict";
    var n = function(t) {
            this.defaultColor = "#ffffff", this.type = "static", this.map = t, this.getColor = function(t) {
                return void 0 !== this.map[t] ? this.map[t] : this.defaultColor
            }
        },
        r = function(t, e) {
            this.type = "dyn", this.opt = e, void 0 !== t.init ? (t.init.call(this), this.getColor = t.run, this.reset = t.init) : this.getColor = t
        };
    t.exports.stat = n, t.exports.dyn = r
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#5858a7",
        R: "#6b6b94",
        N: "#64649b",
        D: "#2121de",
        C: "#9d9d62",
        Q: "#8c8c73",
        E: "#0000ff",
        G: "#4949b6",
        H: "#60609f",
        I: "#ecec13",
        L: "#b2b24d",
        K: "#4747b8",
        M: "#82827d",
        F: "#c2c23d",
        P: "#2323dc",
        S: "#4949b6",
        T: "#9d9d62",
        W: "#c0c03f",
        Y: "#d3d32c",
        V: "#ffff00",
        B: "#4343bc",
        X: "#797986",
        Z: "#4747b8"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#ccff00",
        R: "#0000ff",
        N: "#cc00ff",
        D: "#ff0000",
        C: "#ffff00",
        Q: "#ff00cc",
        E: "#ff0066",
        G: "#ff9900",
        H: "#0066ff",
        I: "#66ff00",
        L: "#33ff00",
        K: "#6600ff",
        M: "#00ff00",
        F: "#00ff66",
        P: "#ffcc00",
        S: "#ff3300",
        T: "#ff6600",
        W: "#00ccff",
        Y: "#00ffcc",
        V: "#99ff00",
        B: "#fff",
        X: "#fff",
        Z: "#fff"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#2cd3d3",
        R: "#708f8f",
        N: "#ff0000",
        D: "#e81717",
        C: "#a85757",
        Q: "#3fc0c0",
        E: "#778888",
        G: "#ff0000",
        H: "#708f8f",
        I: "#00ffff",
        L: "#1ce3e3",
        K: "#7e8181",
        M: "#1ee1e1",
        F: "#1ee1e1",
        P: "#f60909",
        S: "#e11e1e",
        T: "#738c8c",
        W: "#738c8c",
        Y: "#9d6262",
        V: "#07f8f8",
        B: "#f30c0c",
        X: "#7c8383",
        Z: "#5ba4a4"
    }
}, function(t, e) {
    "use strict";
    t.exports = {
        A: "#ffafaf",
        R: "#6464ff",
        N: "#00ff00",
        D: "#ff0000",
        C: "#ffff00",
        Q: "#00ff00",
        E: "#ff0000",
        G: "#ff00ff",
        H: "#6464ff",
        I: "#ffafaf",
        L: "#ffafaf",
        K: "#6464ff",
        M: "#ffafaf",
        F: "#ffc800",
        P: "#ff00ff",
        S: "#00ff00",
        T: "#00ff00",
        W: "#ffc800",
        Y: "#ffc800",
        V: "#ffafaf",
        B: "#fff",
        X: "#fff",
        Z: "#fff"
    }
}, function(t, e, n) {
    var r, i, o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
        return typeof t
    } : function(t) {
        return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
    };
    (function() {
        function n(t) {
            function e(e, n, r, i, o, s) {
                for (; o >= 0 && o < s; o += t) {
                    var u = i ? i[o] : o;
                    r = n(r, e[u], u, e)
                }
                return r
            }
            return function(n, r, i, o) {
                r = k(r, o, 4);
                var s = !C(n) && S.keys(n),
                    u = (s || n).length,
                    a = t > 0 ? 0 : u - 1;
                return arguments.length < 3 && (i = n[s ? s[a] : a], a += t), e(n, r, i, s, a, u)
            }
        }

        function s(t) {
            return function(e, n, r) {
                n = j(n, r);
                for (var i = A(e), o = t > 0 ? 0 : i - 1; o >= 0 && o < i; o += t)
                    if (n(e[o], o, e)) return o;
                return -1
            }
        }

        function u(t, e, n) {
            return function(r, i, o) {
                var s = 0,
                    u = A(r);
                if ("number" == typeof o) t > 0 ? s = o >= 0 ? o : Math.max(o + u, s) : u = o >= 0 ? Math.min(o + 1, u) : o + u + 1;
                else if (n && o && u) return o = n(r, i), r[o] === i ? o : -1;
                if (i !== i) return o = e(g.call(r, s, u), S.isNaN), o >= 0 ? o + s : -1;
                for (o = t > 0 ? s : u - 1; o >= 0 && o < u; o += t)
                    if (r[o] === i) return o;
                return -1
            }
        }

        function a(t, e) {
            var n = R.length,
                r = t.constructor,
                i = S.isFunction(r) && r.prototype || h,
                o = "constructor";
            for (S.has(t, o) && !S.contains(e, o) && e.push(o); n--;) o = R[n], o in t && t[o] !== i[o] && !S.contains(e, o) && e.push(o)
        }
        var l = this,
            c = l._,
            f = Array.prototype,
            h = Object.prototype,
            d = Function.prototype,
            p = f.push,
            g = f.slice,
            v = h.toString,
            m = h.hasOwnProperty,
            y = Array.isArray,
            _ = Object.keys,
            b = d.bind,
            x = Object.create,
            w = function() {},
            S = function G(t) {
                return t instanceof G ? t : this instanceof G ? void(this._wrapped = t) : new G(t)
            };
        "undefined" != typeof t && t.exports && (e = t.exports = S), e._ = S, S.VERSION = "1.8.3";
        var k = function(t, e, n) {
                if (void 0 === e) return t;
                switch (null == n ? 3 : n) {
                    case 1:
                        return function(n) {
                            return t.call(e, n)
                        };
                    case 2:
                        return function(n, r) {
                            return t.call(e, n, r)
                        };
                    case 3:
                        return function(n, r, i) {
                            return t.call(e, n, r, i)
                        };
                    case 4:
                        return function(n, r, i, o) {
                            return t.call(e, n, r, i, o)
                        }
                }
                return function() {
                    return t.apply(e, arguments)
                }
            },
            j = function(t, e, n) {
                return null == t ? S.identity : S.isFunction(t) ? k(t, e, n) : S.isObject(t) ? S.matcher(t) : S.property(t)
            };
        S.iteratee = function(t, e) {
            return j(t, e, 1 / 0)
        };
        var O = function(t, e) {
                return function(n) {
                    var r = arguments.length;
                    if (r < 2 || null == n) return n;
                    for (var i = 1; i < r; i++)
                        for (var o = arguments[i], s = t(o), u = s.length, a = 0; a < u; a++) {
                            var l = s[a];
                            e && void 0 !== n[l] || (n[l] = o[l])
                        }
                    return n
                }
            },
            E = function(t) {
                if (!S.isObject(t)) return {};
                if (x) return x(t);
                w.prototype = t;
                var e = new w;
                return w.prototype = null, e
            },
            M = function(t) {
                return function(e) {
                    return null == e ? void 0 : e[t]
                }
            },
            z = Math.pow(2, 53) - 1,
            A = M("length"),
            C = function(t) {
                var e = A(t);
                return "number" == typeof e && e >= 0 && e <= z
            };
        S.each = S.forEach = function(t, e, n) {
            e = k(e, n);
            var r, i;
            if (C(t))
                for (r = 0, i = t.length; r < i; r++) e(t[r], r, t);
            else {
                var o = S.keys(t);
                for (r = 0, i = o.length; r < i; r++) e(t[o[r]], o[r], t)
            }
            return t
        }, S.map = S.collect = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = Array(i), s = 0; s < i; s++) {
                var u = r ? r[s] : s;
                o[s] = e(t[u], u, t)
            }
            return o
        }, S.reduce = S.foldl = S.inject = n(1), S.reduceRight = S.foldr = n(-1), S.find = S.detect = function(t, e, n) {
            var r;
            if (r = C(t) ? S.findIndex(t, e, n) : S.findKey(t, e, n), void 0 !== r && r !== -1) return t[r]
        }, S.filter = S.select = function(t, e, n) {
            var r = [];
            return e = j(e, n), S.each(t, function(t, n, i) {
                e(t, n, i) && r.push(t)
            }), r
        }, S.reject = function(t, e, n) {
            return S.filter(t, S.negate(j(e)), n)
        }, S.every = S.all = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (!e(t[s], s, t)) return !1
            }
            return !0
        }, S.some = S.any = function(t, e, n) {
            e = j(e, n);
            for (var r = !C(t) && S.keys(t), i = (r || t).length, o = 0; o < i; o++) {
                var s = r ? r[o] : o;
                if (e(t[s], s, t)) return !0
            }
            return !1
        }, S.contains = S.includes = S.include = function(t, e, n, r) {
            return C(t) || (t = S.values(t)), ("number" != typeof n || r) && (n = 0), S.indexOf(t, e, n) >= 0
        }, S.invoke = function(t, e) {
            var n = g.call(arguments, 2),
                r = S.isFunction(e);
            return S.map(t, function(t) {
                var i = r ? e : t[e];
                return null == i ? i : i.apply(t, n)
            })
        }, S.pluck = function(t, e) {
            return S.map(t, S.property(e))
        }, S.where = function(t, e) {
            return S.filter(t, S.matcher(e))
        }, S.findWhere = function(t, e) {
            return S.find(t, S.matcher(e))
        }, S.max = function(t, e, n) {
            var r, i, o = -(1 / 0),
                s = -(1 / 0);
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r > o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i > s || i === -(1 / 0) && o === -(1 / 0)) && (o = t, s = i)
            });
            return o
        }, S.min = function(t, e, n) {
            var r, i, o = 1 / 0,
                s = 1 / 0;
            if (null == e && null != t) {
                t = C(t) ? t : S.values(t);
                for (var u = 0, a = t.length; u < a; u++) r = t[u], r < o && (o = r)
            } else e = j(e, n), S.each(t, function(t, n, r) {
                i = e(t, n, r), (i < s || i === 1 / 0 && o === 1 / 0) && (o = t, s = i)
            });
            return o
        }, S.shuffle = function(t) {
            for (var e, n = C(t) ? t : S.values(t), r = n.length, i = Array(r), o = 0; o < r; o++) e = S.random(0, o), e !== o && (i[o] = i[e]), i[e] = n[o];
            return i
        }, S.sample = function(t, e, n) {
            return null == e || n ? (C(t) || (t = S.values(t)), t[S.random(t.length - 1)]) : S.shuffle(t).slice(0, Math.max(0, e))
        }, S.sortBy = function(t, e, n) {
            return e = j(e, n), S.pluck(S.map(t, function(t, n, r) {
                return {
                    value: t,
                    index: n,
                    criteria: e(t, n, r)
                }
            }).sort(function(t, e) {
                var n = t.criteria,
                    r = e.criteria;
                if (n !== r) {
                    if (n > r || void 0 === n) return 1;
                    if (n < r || void 0 === r) return -1
                }
                return t.index - e.index
            }), "value")
        };
        var T = function(t) {
            return function(e, n, r) {
                var i = {};
                return n = j(n, r), S.each(e, function(r, o) {
                    t(i, r, n(r, o, e))
                }), i
            }
        };
        S.groupBy = T(function(t, e, n) {
            S.has(t, n) ? t[n].push(e) : t[n] = [e]
        }), S.indexBy = T(function(t, e, n) {
            t[n] = e
        }), S.countBy = T(function(t, e, n) {
            S.has(t, n) ? t[n]++ : t[n] = 1
        }), S.toArray = function(t) {
            return t ? S.isArray(t) ? g.call(t) : C(t) ? S.map(t, S.identity) : S.values(t) : []
        }, S.size = function(t) {
            return null == t ? 0 : C(t) ? t.length : S.keys(t).length
        }, S.partition = function(t, e, n) {
            e = j(e, n);
            var r = [],
                i = [];
            return S.each(t, function(t, n, o) {
                (e(t, n, o) ? r : i).push(t)
            }), [r, i]
        }, S.first = S.head = S.take = function(t, e, n) {
            if (null != t) return null == e || n ? t[0] : S.initial(t, t.length - e)
        }, S.initial = function(t, e, n) {
            return g.call(t, 0, Math.max(0, t.length - (null == e || n ? 1 : e)))
        }, S.last = function(t, e, n) {
            if (null != t) return null == e || n ? t[t.length - 1] : S.rest(t, Math.max(0, t.length - e))
        }, S.rest = S.tail = S.drop = function(t, e, n) {
            return g.call(t, null == e || n ? 1 : e)
        }, S.compact = function(t) {
            return S.filter(t, S.identity)
        };
        var I = function K(t, e, n, r) {
            for (var i = [], o = 0, s = r || 0, u = A(t); s < u; s++) {
                var a = t[s];
                if (C(a) && (S.isArray(a) || S.isArguments(a))) {
                    e || (a = K(a, e, n));
                    var l = 0,
                        c = a.length;
                    for (i.length += c; l < c;) i[o++] = a[l++]
                } else n || (i[o++] = a)
            }
            return i
        };
        S.flatten = function(t, e) {
            return I(t, e, !1)
        }, S.without = function(t) {
            return S.difference(t, g.call(arguments, 1))
        }, S.uniq = S.unique = function(t, e, n, r) {
            S.isBoolean(e) || (r = n, n = e, e = !1), null != n && (n = j(n, r));
            for (var i = [], o = [], s = 0, u = A(t); s < u; s++) {
                var a = t[s],
                    l = n ? n(a, s, t) : a;
                e ? (s && o === l || i.push(a), o = l) : n ? S.contains(o, l) || (o.push(l), i.push(a)) : S.contains(i, a) || i.push(a)
            }
            return i
        }, S.union = function() {
            return S.uniq(I(arguments, !0, !0))
        }, S.intersection = function(t) {
            for (var e = [], n = arguments.length, r = 0, i = A(t); r < i; r++) {
                var o = t[r];
                if (!S.contains(e, o)) {
                    for (var s = 1; s < n && S.contains(arguments[s], o); s++);
                    s === n && e.push(o)
                }
            }
            return e
        }, S.difference = function(t) {
            var e = I(arguments, !0, !0, 1);
            return S.filter(t, function(t) {
                return !S.contains(e, t)
            })
        }, S.zip = function() {
            return S.unzip(arguments)
        }, S.unzip = function(t) {
            for (var e = t && S.max(t, A).length || 0, n = Array(e), r = 0; r < e; r++) n[r] = S.pluck(t, r);
            return n
        }, S.object = function(t, e) {
            for (var n = {}, r = 0, i = A(t); r < i; r++) e ? n[t[r]] = e[r] : n[t[r][0]] = t[r][1];
            return n
        }, S.findIndex = s(1), S.findLastIndex = s(-1), S.sortedIndex = function(t, e, n, r) {
            n = j(n, r, 1);
            for (var i = n(e), o = 0, s = A(t); o < s;) {
                var u = Math.floor((o + s) / 2);
                n(t[u]) < i ? o = u + 1 : s = u
            }
            return o
        }, S.indexOf = u(1, S.findIndex, S.sortedIndex), S.lastIndexOf = u(-1, S.findLastIndex), S.range = function(t, e, n) {
            null == e && (e = t || 0, t = 0), n = n || 1;
            for (var r = Math.max(Math.ceil((e - t) / n), 0), i = Array(r), o = 0; o < r; o++, t += n) i[o] = t;
            return i
        };
        var N = function(t, e, n, r, i) {
            if (!(r instanceof e)) return t.apply(n, i);
            var o = E(t.prototype),
                s = t.apply(o, i);
            return S.isObject(s) ? s : o
        };
        S.bind = function(t, e) {
            if (b && t.bind === b) return b.apply(t, g.call(arguments, 1));
            if (!S.isFunction(t)) throw new TypeError("Bind must be called on a function");
            var n = g.call(arguments, 2);
            return function r() {
                return N(t, r, e, this, n.concat(g.call(arguments)))
            }
        }, S.partial = function(t) {
            var e = g.call(arguments, 1);
            return function n() {
                for (var r = 0, i = e.length, o = Array(i), s = 0; s < i; s++) o[s] = e[s] === S ? arguments[r++] : e[s];
                for (; r < arguments.length;) o.push(arguments[r++]);
                return N(t, n, this, this, o)
            }
        }, S.bindAll = function(t) {
            var e, n, r = arguments.length;
            if (r <= 1) throw new Error("bindAll must be passed function names");
            for (e = 1; e < r; e++) n = arguments[e], t[n] = S.bind(t[n], t);
            return t
        }, S.memoize = function(t, e) {
            var n = function r(n) {
                var i = r.cache,
                    o = "" + (e ? e.apply(this, arguments) : n);
                return S.has(i, o) || (i[o] = t.apply(this, arguments)), i[o]
            };
            return n.cache = {}, n
        }, S.delay = function(t, e) {
            var n = g.call(arguments, 2);
            return setTimeout(function() {
                return t.apply(null, n)
            }, e)
        }, S.defer = S.partial(S.delay, S, 1), S.throttle = function(t, e, n) {
            var r, i, o, s = null,
                u = 0;
            n || (n = {});
            var a = function() {
                u = n.leading === !1 ? 0 : S.now(), s = null, o = t.apply(r, i), s || (r = i = null)
            };
            return function() {
                var l = S.now();
                u || n.leading !== !1 || (u = l);
                var c = e - (l - u);
                return r = this, i = arguments, c <= 0 || c > e ? (s && (clearTimeout(s), s = null), u = l, o = t.apply(r, i), s || (r = i = null)) : s || n.trailing === !1 || (s = setTimeout(a, c)), o
            }
        }, S.debounce = function(t, e, n) {
            var r, i, o, s, u, a = function l() {
                var a = S.now() - s;
                a < e && a >= 0 ? r = setTimeout(l, e - a) : (r = null, n || (u = t.apply(o, i), r || (o = i = null)))
            };
            return function() {
                o = this, i = arguments, s = S.now();
                var l = n && !r;
                return r || (r = setTimeout(a, e)), l && (u = t.apply(o, i), o = i = null), u
            }
        }, S.wrap = function(t, e) {
            return S.partial(e, t)
        }, S.negate = function(t) {
            return function() {
                return !t.apply(this, arguments)
            }
        }, S.compose = function() {
            var t = arguments,
                e = t.length - 1;
            return function() {
                for (var n = e, r = t[e].apply(this, arguments); n--;) r = t[n].call(this, r);
                return r
            }
        }, S.after = function(t, e) {
            return function() {
                if (--t < 1) return e.apply(this, arguments)
            }
        }, S.before = function(t, e) {
            var n;
            return function() {
                return --t > 0 && (n = e.apply(this, arguments)), t <= 1 && (e = null), n
            }
        }, S.once = S.partial(S.before, 2);
        var L = !{
                toString: null
            }.propertyIsEnumerable("toString"),
            R = ["valueOf", "isPrototypeOf", "toString", "propertyIsEnumerable", "hasOwnProperty", "toLocaleString"];
        S.keys = function(t) {
            if (!S.isObject(t)) return [];
            if (_) return _(t);
            var e = [];
            for (var n in t) S.has(t, n) && e.push(n);
            return L && a(t, e), e
        }, S.allKeys = function(t) {
            if (!S.isObject(t)) return [];
            var e = [];
            for (var n in t) e.push(n);
            return L && a(t, e), e
        }, S.values = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = t[e[i]];
            return r
        }, S.mapObject = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = i.length, s = {}, u = 0; u < o; u++) r = i[u], s[r] = e(t[r], r, t);
            return s
        }, S.pairs = function(t) {
            for (var e = S.keys(t), n = e.length, r = Array(n), i = 0; i < n; i++) r[i] = [e[i], t[e[i]]];
            return r
        }, S.invert = function(t) {
            for (var e = {}, n = S.keys(t), r = 0, i = n.length; r < i; r++) e[t[n[r]]] = n[r];
            return e
        }, S.functions = S.methods = function(t) {
            var e = [];
            for (var n in t) S.isFunction(t[n]) && e.push(n);
            return e.sort()
        }, S.extend = O(S.allKeys), S.extendOwn = S.assign = O(S.keys), S.findKey = function(t, e, n) {
            e = j(e, n);
            for (var r, i = S.keys(t), o = 0, s = i.length; o < s; o++)
                if (r = i[o], e(t[r], r, t)) return r
        }, S.pick = function(t, e, n) {
            var r, i, o = {},
                s = t;
            if (null == s) return o;
            S.isFunction(e) ? (i = S.allKeys(s), r = k(e, n)) : (i = I(arguments, !1, !1, 1), r = function(t, e, n) {
                return e in n
            }, s = Object(s));
            for (var u = 0, a = i.length; u < a; u++) {
                var l = i[u],
                    c = s[l];
                r(c, l, s) && (o[l] = c)
            }
            return o
        }, S.omit = function(t, e, n) {
            if (S.isFunction(e)) e = S.negate(e);
            else {
                var r = S.map(I(arguments, !1, !1, 1), String);
                e = function(t, e) {
                    return !S.contains(r, e)
                }
            }
            return S.pick(t, e, n)
        }, S.defaults = O(S.allKeys, !0), S.create = function(t, e) {
            var n = E(t);
            return e && S.extendOwn(n, e), n
        }, S.clone = function(t) {
            return S.isObject(t) ? S.isArray(t) ? t.slice() : S.extend({}, t) : t
        }, S.tap = function(t, e) {
            return e(t), t
        }, S.isMatch = function(t, e) {
            var n = S.keys(e),
                r = n.length;
            if (null == t) return !r;
            for (var i = Object(t), o = 0; o < r; o++) {
                var s = n[o];
                if (e[s] !== i[s] || !(s in i)) return !1
            }
            return !0
        };
        var q = function X(t, e, n, r) {
            if (t === e) return 0 !== t || 1 / t === 1 / e;
            if (null == t || null == e) return t === e;
            t instanceof S && (t = t._wrapped), e instanceof S && (e = e._wrapped);
            var i = v.call(t);
            if (i !== v.call(e)) return !1;
            switch (i) {
                case "[object RegExp]":
                case "[object String]":
                    return "" + t == "" + e;
                case "[object Number]":
                    return +t !== +t ? +e !== +e : 0 === +t ? 1 / +t === 1 / e : +t === +e;
                case "[object Date]":
                case "[object Boolean]":
                    return +t === +e
            }
            var s = "[object Array]" === i;
            if (!s) {
                if ("object" != ("undefined" == typeof t ? "undefined" : o(t)) || "object" != ("undefined" == typeof e ? "undefined" : o(e))) return !1;
                var u = t.constructor,
                    a = e.constructor;
                if (u !== a && !(S.isFunction(u) && u instanceof u && S.isFunction(a) && a instanceof a) && "constructor" in t && "constructor" in e) return !1
            }
            n = n || [], r = r || [];
            for (var l = n.length; l--;)
                if (n[l] === t) return r[l] === e;
            if (n.push(t), r.push(e), s) {
                if (l = t.length, l !== e.length) return !1;
                for (; l--;)
                    if (!X(t[l], e[l], n, r)) return !1
            } else {
                var c, f = S.keys(t);
                if (l = f.length, S.keys(e).length !== l) return !1;
                for (; l--;)
                    if (c = f[l], !S.has(e, c) || !X(t[c], e[c], n, r)) return !1
            }
            return n.pop(), r.pop(), !0
        };
        S.isEqual = function(t, e) {
            return q(t, e)
        }, S.isEmpty = function(t) {
            return null == t || (C(t) && (S.isArray(t) || S.isString(t) || S.isArguments(t)) ? 0 === t.length : 0 === S.keys(t).length)
        }, S.isElement = function(t) {
            return !(!t || 1 !== t.nodeType)
        }, S.isArray = y || function(t) {
            return "[object Array]" === v.call(t)
        }, S.isObject = function(t) {
            var e = "undefined" == typeof t ? "undefined" : o(t);
            return "function" === e || "object" === e && !!t
        }, S.each(["Arguments", "Function", "String", "Number", "Date", "RegExp", "Error"], function(t) {
            S["is" + t] = function(e) {
                return v.call(e) === "[object " + t + "]"
            }
        }), S.isArguments(arguments) || (S.isArguments = function(t) {
            return S.has(t, "callee")
        }), "function" != typeof /./ && "object" != ("undefined" == typeof Int8Array ? "undefined" : o(Int8Array)) && (S.isFunction = function(t) {
            return "function" == typeof t || !1
        }), S.isFinite = function(t) {
            return isFinite(t) && !isNaN(parseFloat(t))
        }, S.isNaN = function(t) {
            return S.isNumber(t) && t !== +t
        }, S.isBoolean = function(t) {
            return t === !0 || t === !1 || "[object Boolean]" === v.call(t)
        }, S.isNull = function(t) {
            return null === t
        }, S.isUndefined = function(t) {
            return void 0 === t
        }, S.has = function(t, e) {
            return null != t && m.call(t, e)
        }, S.noConflict = function() {
            return l._ = c, this
        }, S.identity = function(t) {
            return t
        }, S.constant = function(t) {
            return function() {
                return t
            }
        }, S.noop = function() {}, S.property = M, S.propertyOf = function(t) {
            return null == t ? function() {} : function(e) {
                return t[e]
            }
        }, S.matcher = S.matches = function(t) {
            return t = S.extendOwn({}, t),
                function(e) {
                    return S.isMatch(e, t)
                }
        }, S.times = function(t, e, n) {
            var r = Array(Math.max(0, t));
            e = k(e, n, 1);
            for (var i = 0; i < t; i++) r[i] = e(i);
            return r
        }, S.random = function(t, e) {
            return null == e && (e = t, t = 0), t + Math.floor(Math.random() * (e - t + 1))
        }, S.now = Date.now || function() {
            return (new Date).getTime()
        };
        var F = {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': "&quot;",
                "'": "&#x27;",
                "`": "&#x60;"
            },
            P = S.invert(F),
            B = function(t) {
                var e = function(e) {
                        return t[e]
                    },
                    n = "(?:" + S.keys(t).join("|") + ")",
                    r = RegExp(n),
                    i = RegExp(n, "g");
                return function(t) {
                    return t = null == t ? "" : "" + t, r.test(t) ? t.replace(i, e) : t
                }
            };
        S.escape = B(F), S.unescape = B(P), S.result = function(t, e, n) {
            var r = null == t ? void 0 : t[e];
            return void 0 === r && (r = n), S.isFunction(r) ? r.call(t) : r
        };
        var W = 0;
        S.uniqueId = function(t) {
            var e = ++W + "";
            return t ? t + e : e
        }, S.templateSettings = {
            evaluate: /<%([\s\S]+?)%>/g,
            interpolate: /<%=([\s\S]+?)%>/g,
            escape: /<%-([\s\S]+?)%>/g
        };
        var D = /(.)^/,
            H = {
                "'": "'",
                "\\": "\\",
                "\r": "r",
                "\n": "n",
                "\u2028": "u2028",
                "\u2029": "u2029"
            },
            U = /\\|'|\r|\n|\u2028|\u2029/g,
            V = function(t) {
                return "\\" + H[t]
            };
        S.template = function(t, e, n) {
            !e && n && (e = n), e = S.defaults({}, e, S.templateSettings);
            var r = RegExp([(e.escape || D).source, (e.interpolate || D).source, (e.evaluate || D).source].join("|") + "|$", "g"),
                i = 0,
                o = "__p+='";
            t.replace(r, function(e, n, r, s, u) {
                return o += t.slice(i, u).replace(U, V), i = u + e.length, n ? o += "'+\n((__t=(" + n + "))==null?'':_.escape(__t))+\n'" : r ? o += "'+\n((__t=(" + r + "))==null?'':__t)+\n'" : s && (o += "';\n" + s + "\n__p+='"), e
            }), o += "';\n", e.variable || (o = "with(obj||{}){\n" + o + "}\n"), o = "var __t,__p='',__j=Array.prototype.join,print=function(){__p+=__j.call(arguments,'');};\n" + o + "return __p;\n";
            try {
                var s = new Function(e.variable || "obj", "_", o)
            } catch (u) {
                throw u.source = o, u
            }
            var a = function(t) {
                    return s.call(this, t, S)
                },
                l = e.variable || "obj";
            return a.source = "function(" + l + "){\n" + o + "}", a
        }, S.chain = function(t) {
            var e = S(t);
            return e._chain = !0, e
        };
        var $ = function(t, e) {
            return t._chain ? S(e).chain() : e
        };
        S.mixin = function(t) {
            S.each(S.functions(t), function(e) {
                var n = S[e] = t[e];
                S.prototype[e] = function() {
                    var t = [this._wrapped];
                    return p.apply(t, arguments), $(this, n.apply(S, t))
                }
            })
        }, S.mixin(S), S.each(["pop", "push", "reverse", "shift", "sort", "splice", "unshift"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                var n = this._wrapped;
                return e.apply(n, arguments), "shift" !== t && "splice" !== t || 0 !== n.length || delete n[0], $(this, n)
            }
        }), S.each(["concat", "join", "slice"], function(t) {
            var e = f[t];
            S.prototype[t] = function() {
                return $(this, e.apply(this._wrapped, arguments))
            }
        }), S.prototype.value = function() {
            return this._wrapped
        }, S.prototype.valueOf = S.prototype.toJSON = S.prototype.value, S.prototype.toString = function() {
            return "" + this._wrapped
        }, r = [], i = function() {
            return S
        }.apply(e, r), !(void 0 !== i && (t.exports = i))
    }).call(void 0)
}, function(t, e, n) {
    "use strict";

    function r(t) {
        if (t && t.__esModule) return t;
        var e = {};
        if (null != t)
            for (var n in t) Object.prototype.hasOwnProperty.call(t, n) && (e[n] = t[n]);
        return e["default"] = t, e
    }

    function i(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    }), e.version = e.io = e.$ = e.boneView = e.view = e.selcol = e.selection = e.utils = e.menu = e.model = e.msa = void 0;
    var o = n(8);
    Object.defineProperty(e, "selection", {
        enumerable: !0,
        get: function() {
            return i(o)["default"]
        }
    });
    var s = n(9);
    Object.defineProperty(e, "selcol", {
        enumerable: !0,
        get: function() {
            return i(s)["default"]
        }
    });
    var u = n(2);
    Object.defineProperty(e, "view", {
        enumerable: !0,
        get: function() {
            return i(u)["default"]
        }
    });
    var a = n(4);
    Object.defineProperty(e, "boneView", {
        enumerable: !0,
        get: function() {
            return i(a)["default"]
        }
    });
    var l = n(5);
    Object.defineProperty(e, "$", {
        enumerable: !0,
        get: function() {
            return i(l)["default"]
        }
    });
    var c = n(125),
        f = i(c),
        h = n(51),
        d = r(h),
        p = n(113),
        g = r(p),
        v = n(37),
        m = r(v),
        y = n(10),
        _ = function() {
            var t = function(t) {
                return f["default"].apply(this, t)
            };
            return t.prototype = f["default"].prototype, new t(arguments)
        };
    e["default"] = _, e.msa = f["default"], e.model = d, e.menu = g, e.utils = m;
    var b = {
        xhr: n(21),
        fasta: y.fasta,
        clustal: y.clustal,
        gff: y.gff
    };
    e.io = b;
    var x = "imported";
    x = "1.0.3", e.version = x
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(121),
        o = r(i),
        s = n(119),
        u = r(s),
        a = n(123),
        l = r(a),
        c = n(124),
        f = r(c),
        h = n(115),
        d = r(h),
        p = n(122),
        g = r(p),
        v = n(118),
        m = r(v),
        y = n(117),
        _ = r(y),
        b = n(120),
        x = r(b),
        w = n(116),
        S = r(w),
        k = n(114),
        j = r(k),
        O = n(4),
        E = O.extend({
            initialize: function(t) {
                if (!t.msa) throw new Error("No msa instance provided. Please provide .msa");
                if (this.msa = t.msa, this.msa.g.menuconfig = new j["default"](t.menu), this.addView("10_import", new o["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g,
                        msa: this.msa
                    })), this.addView("15_ordering", new g["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g
                    })), this.addView("20_filter", new u["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g
                    })), this.addView("30_selection", new l["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g
                    })), this.addView("40_vis", new f["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g
                    })), this.addView("50_color", new d["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g
                    })), this.addView("70_extra", new m["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g,
                        msa: this.msa
                    })), this.addView("80_export", new _["default"]({
                        model: this.msa.seqs,
                        g: this.msa.g,
                        msa: this.msa
                    })), this.addView("90_help", new x["default"]({
                        g: this.msa.g
                    })), this.msa.g.config.get("debug")) return this.addView("95_debug", new S["default"]({
                    g: this.msa.g
                }))
            },
            render: function() {
                return this.renderSubviews(), this.el.setAttribute("class", "smenubar"), this.el.appendChild(document.createElement("p"))
            }
        });
    e["default"] = E
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(112);
    Object.defineProperty(e, "defaultmenu", {
        enumerable: !0,
        get: function() {
            return r(i)["default"]
        }
    });
    var o = n(6);
    Object.defineProperty(e, "menubuilder", {
        enumerable: !0,
        get: function() {
            return r(o)["default"]
        }
    })
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(1).Model,
        i = r.extend({
            constructor: function(t, e) {
                return "small" == t && (t = this.small), r.apply(this, [t])
            },
            small: {
                menuFontsize: "12px"
            },
            defaults: {
                menuFontsize: "14px",
                menuItemFontsize: "14px",
                menuItemLineHeight: "14px",
                menuMarginLeft: "3px",
                menuPadding: "3px 4px 3px 4px"
            }
        });
    e["default"] = i
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = n(7),
        u = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block",
                    this.listenTo(this.g.colorscheme, "change", function() {
                        return this.render()
                    })
            },
            render: function() {
                var t = this.setName("Color scheme");
                this.removeAllNodes();
                for (var e, n = this.getColorschemes(), r = 0; r < n.length; r++) e = n[r], this.addScheme(t, e);
                return this.grey(t), s.removeAllChilds(this.el), this.el.appendChild(this.buildDOM()), this
            },
            addScheme: function(t, e) {
                var n = this,
                    r = {};
                return this.g.colorscheme.get("scheme") === e.id && (r.backgroundColor = "#77ED80"), this.addNode(e.name, function() {
                    n.g.colorscheme.set("scheme", e.id)
                }, {
                    style: r
                })
            },
            getColorschemes: function() {
                var t = [];
                return t.push({
                    name: "Taylor",
                    id: "taylor"
                }), t.push({
                    name: "Buried",
                    id: "buried"
                }), t.push({
                    name: "Cinema",
                    id: "cinema"
                }), t.push({
                    name: "Clustal",
                    id: "clustal"
                }), t.push({
                    name: "Clustal2",
                    id: "clustal2"
                }), t.push({
                    name: "Helix",
                    id: "helix"
                }), t.push({
                    name: "Hydrophobicity",
                    id: "hydro"
                }), t.push({
                    name: "Lesk",
                    id: "lesk"
                }), t.push({
                    name: "MAE",
                    id: "mae"
                }), t.push({
                    name: "Nucleotide",
                    id: "nucleotide"
                }), t.push({
                    name: "Purine",
                    id: "purine"
                }), t.push({
                    name: "PID",
                    id: "pid"
                }), t.push({
                    name: "Strand",
                    id: "strand"
                }), t.push({
                    name: "Turn",
                    id: "turn"
                }), t.push({
                    name: "Zappo",
                    id: "zappo"
                }), t.push({
                    name: "No color",
                    id: "foo"
                }), t
            },
            grey: function(t) {}
        });
    e["default"] = u
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block"
            },
            render: function() {
                var t = this;
                return this.setName("Debug"), this.addNode("Get the code", function() {
                    return window.open("https://github.com/wilzbach/msa")
                }), this.addNode("Toggle mouseover events", function() {
                    return t.g.config.set("registerMouseHover", !t.g.config.get("registerMouseHover")), t.g.onAll(function() {
                        return console.log(arguments)
                    })
                }), this.addNode("Minimized width", function() {
                    return t.g.zoomer.set("alignmentWidth", 600)
                }), this.addNode("Minimized height", function() {
                    return t.g.zoomer.set("alignmentHeight", 120)
                }), this.el.appendChild(this.buildDOM()), this
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = n(16),
        u = r(s),
        a = n(10),
        l = (a.fasta.write, "↪"),
        c = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.msa = t.msa, this.el.style.display = "inline-block"
            },
            render: function() {
                var t = this;
                return this.setName("Export"), this.addNode("Share view (URL)" + l, function() {
                    return u["default"].shareLink(t.msa, function(t) {
                        return window.open(t, "_blank")
                    })
                }), this.addNode("View in Jalview", function() {
                    var e = t.g.config.get("url");
                    return "undefined" == typeof e || null === e ? alert("Sequence weren't imported via an URL") : e.indexOf("localhost") ? u["default"].publishWeb(t.msa, function(e) {
                        return u["default"].openInJalview(e, t.g.colorscheme.get("scheme"))
                    }) : u["default"].openInJalview(e, t.g.colorscheme.get("scheme"))
                }), this.addNode("Export alignment (FASTA)", function() {
                    return u["default"].saveAsFile(t.msa, "all.fasta")
                }), this.addNode("Export alignment (URL)", function() {
                    return u["default"].publishWeb(t.msa, function(t) {
                        return window.open(t, "_blank")
                    })
                }), this.addNode("Export selected sequences (FASTA)", function() {
                    return u["default"].saveSelection(t.msa, "selection.fasta")
                }), this.addNode("Export features (GFF)", function() {
                    return u["default"].saveAnnots(t.msa, "features.gff3")
                }), this.addNode("Export MSA image (PNG)", function() {
                    return u["default"].saveAsImg(t.msa, "biojs-msa.png")
                }), this.el.appendChild(this.buildDOM()), this
            }
        });
    e["default"] = c
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = n(36),
        u = r(s),
        a = n(18),
        l = (r(a), n(21), o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block", this.msa = t.msa
            },
            render: function() {
                var t = this;
                this.setName("Extras");
                var e = this.g.stats;
                return this.msa, this.addNode("Add consensus seq", function() {
                    var n = e.consensus(),
                        r = new u["default"]({
                            seq: n,
                            id: "0c",
                            name: "Consenus"
                        });
                    return t.model.add(r), t.model.setRef(r), t.model.comparator = function(t) {
                        return !t.get("ref")
                    }, t.model.sort()
                }), this.addNode("Jump to a column", function() {
                    var e = prompt("Column", "20");
                    return e < 0 || e > t.model.getMaxLength() || isNaN(e) ? void alert("invalid column") : t.g.zoomer.setLeftOffset(e)
                }), this.el.appendChild(this.buildDOM()), this
            }
        }));
    e["default"] = l
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block"
            },
            render: function() {
                var t = this;
                return this.setName("Filter"), this.addNode("Hide columns by threshold", function(e) {
                    var n = prompt("Enter threshold (in percent)", 20);
                    n /= 100;
                    for (var r = t.model.getMaxLength(), i = [], o = t.g.stats.scale(t.g.stats.conservation()), s = r - 1, u = 0; 0 < s ? u <= s : u >= s; 0 < s ? u++ : u--) o[u] < n && i.push(u);
                    return t.g.columns.set("hidden", i)
                }), this.addNode("Hide columns by selection", function() {
                    var e = t.g.columns.get("hidden"),
                        n = e.concat(t.g.selcol.getAllColumnBlocks({
                            maxLen: t.model.getMaxLength(),
                            withPos: !0
                        }));
                    return t.g.selcol.reset([]), t.g.columns.set("hidden", n)
                }), this.addNode("Hide columns by gaps", function() {
                    var e = prompt("Enter threshold (in percent)", 20);
                    e /= 100;
                    for (var n = t.model.getMaxLength(), r = [], i = n - 1, o = 0; 0 < i ? o <= i : o >= i; 0 < i ? o++ : o--) {
                        var s = 0,
                            u = 0;
                        t.model.each(function(t) {
                            return "-" === t.get("seq")[o] && s++, u++
                        }), s / u > e && r.push(o)
                    }
                    return t.g.columns.set("hidden", r)
                }), this.addNode("Hide seqs by identity", function() {
                    var e = prompt("Enter threshold (in percent)", 20);
                    return e /= 100, t.model.each(function(t) {
                        if (t.get("identity") < e) return t.set("hidden", !0)
                    })
                }), this.addNode("Hide seqs by selection", function() {
                    var e = t.g.selcol.where({
                            type: "row"
                        }),
                        n = e.map(function(t) {
                            return t.get("seqId")
                        });
                    return t.g.selcol.reset([]), t.model.each(function(t) {
                        if (n.indexOf(t.get("id")) >= 0) return t.set("hidden", !0)
                    })
                }), this.addNode("Hide seqs by gaps", function() {
                    var e = prompt("Enter threshold (in percent)", 40);
                    return t.model.each(function(t, n) {
                        if (t.get("seq").reduce(function(t, e) {
                                return "-" === e ? ++t : void 0
                            }, 0) > e) return t.set("hidden", !0)
                    })
                }), this.addNode("Reset", function() {
                    return t.g.columns.set("hidden", []), t.model.each(function(t) {
                        if (t.get("hidden")) return t.set("hidden", !1)
                    })
                }), this.el.appendChild(this.buildDOM()), this
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g
            },
            render: function() {
                return this.setName("Help"), this.addNode("About the project", function() {
                    return window.open("https://github.com/wilzbach/msa")
                }), this.addNode("Report issues", function() {
                    return window.open("https://github.com/wilzbach/msa/issues")
                }), this.addNode("User manual", function() {
                    return window.open("https://github.com/wilzbach/msa/wiki/User-manual")
                }), this.el.style.display = "inline-block", this.el.appendChild(this.buildDOM()), this
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = n(14),
        u = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block", this.msa = t.msa
            },
            render: function() {
                var t = this,
                    e = this.msa,
                    n = s.mk("input");
                n.type = "file", n.style.display = "none", n.multiple = !0, n.addEventListener("change", function() {
                    var t = n.files || [];
                    return e.u.file.importFiles(t)
                }), this.el.appendChild(n);
                var r = "(Fasta, Clustal, GFF, Jalview features, Newick)";
                return this.setName("Import"), this.addNode("URL", function(e) {
                    var n = prompt("URL " + r, "http://rostlab.org/~goldberg/clustalw2-I20140818-215249-0556-53699878-pg.clustalw");
                    if (n.length > 5) return t.msa.u.file.importURL(n, function() {})
                }), this.addNode("From file " + r, function() {
                    return n.click()
                }), this.addNode("Drag & Drop", function() {
                    return alert("Yep. Just drag & drop your file " + r)
                }), this.el.appendChild(this.buildDOM()), this
            }
        });
    e["default"] = u
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = n(7),
        u = "↑",
        a = "↓",
        l = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.order = "ID", this.el.style.display = "inline-block"
            },
            setOrder: function(t) {
                return this.order = t, this.render()
            },
            render: function() {
                this.setName("Sorting"), this.removeAllNodes();
                for (var t, e = this.getComparators(), n = 0; n < e.length; n++) t = e[n], this._addNode(t);
                var r = this.buildDOM();
                return s.removeAllChilds(this.el), this.el.appendChild(r), this
            },
            _addNode: function(t) {
                var e = this,
                    n = t.text,
                    r = {};
                return n === this.order && (r.backgroundColor = "#77ED80"), this.addNode(n, function() {
                    return null != t.precode && t.precode(), e.model.comparator = t.comparator, e.model.sort(), e.setOrder(t.text)
                }, {
                    style: r
                })
            },
            getComparators: function() {
                var t = this,
                    e = [];
                e.push({
                    text: "ID " + u,
                    comparator: "id"
                }), e.push({
                    text: "ID " + a,
                    comparator: function(t, e) {
                        return -("" + t.get("id")).localeCompare("" + e.get("id"), [], {
                            numeric: !0
                        })
                    }
                }), e.push({
                    text: "Label " + u,
                    comparator: "name"
                }), e.push({
                    text: "Label " + a,
                    comparator: function(t, e) {
                        return -t.get("name").localeCompare(e.get("name"))
                    }
                }), e.push({
                    text: "Seq " + u,
                    comparator: "seq"
                }), e.push({
                    text: "Seq " + a,
                    comparator: function(t, e) {
                        return -t.get("seq").localeCompare(e.get("seq"))
                    }
                });
                var n = function() {
                        return t.ident = t.g.stats.identity()
                    },
                    r = function() {
                        return t.gaps = {}, t.model.each(function(e) {
                            var n = e.attributes.seq;
                            return t.gaps[e.id] = (n.reduce(function(t, e) {
                                return "-" === e ? ++t : void 0
                            }), 0 / n.length)
                        })
                    };
                return e.push({
                    text: "Identity " + u,
                    comparator: function(e, n) {
                        var r = t.ident[e.id] - t.ident[n.id];
                        return r > 0 ? 1 : r < 0 ? -1 : 0
                    },
                    precode: n
                }), e.push({
                    text: "Identity " + a,
                    comparator: function(e, n) {
                        var r = t.ident[e.id] - t.ident[n.id];
                        return r > 0 ? -1 : r < 0 ? 1 : 0
                    },
                    precode: n
                }), e.push({
                    text: "Gaps " + u,
                    comparator: function(e, n) {
                        var r = t.gaps[e.id] - t.gaps[n.id];
                        return r > 0 ? 1 : r < 0 ? -1 : 0
                    },
                    precode: r
                }), e.push({
                    text: "Gaps " + a,
                    comparator: function(e, n) {
                        var r = t.gaps[e.id] - t.gaps[n.id];
                        return r < 0 ? 1 : r > 0 ? -1 : 0
                    },
                    precode: r
                }), e.push({
                    text: "Consensus to top",
                    comparator: function(t) {
                        return !t.get("ref")
                    }
                }), e
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block"
            },
            render: function() {
                var t = this;
                return this.setName("Selection"), this.addNode("Find Motif (supports RegEx)", function() {
                    var e = prompt("your search", "D");
                    return t.g.user.set("searchText", e)
                }), this.addNode("Invert columns", function() {
                    return t.g.selcol.invertCol(function() {
                        var e = [],
                            n = t.model.getMaxLength(),
                            r = 0;
                        if (0 <= n)
                            for (; r <= n;) e.push(r++);
                        else
                            for (; r >= n;) e.push(r--);
                        return e
                    }())
                }), this.addNode("Invert rows", function() {
                    return t.g.selcol.invertRow(t.model.pluck("id"))
                }), this.addNode("Reset", function() {
                    return t.g.selcol.reset()
                }), this.el.appendChild(this.buildDOM()), this
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(6),
        o = r(i),
        s = n(7),
        u = o["default"].extend({
            initialize: function(t) {
                return this.g = t.g, this.el.style.display = "inline-block", this.listenTo(this.g.vis, "change", this.render)
            },
            render: function() {
                var t = this;
                this.removeAllNodes(), this.setName("Vis.elements");
                for (var e, n = this.getVisElements(), r = 0; r < n.length; r++) e = n[r], this._addVisEl(e);
                return this.addNode("Reset", function() {
                    return t.g.vis.set("labels", !0), t.g.vis.set("sequences", !0), t.g.vis.set("metacell", !0), t.g.vis.set("conserv", !0), t.g.vis.set("labelId", !0), t.g.vis.set("labelName", !0), t.g.vis.set("labelCheckbox", !1), t.g.vis.set("seqlogo", !1), t.g.vis.set("gapHeader", !1), t.g.vis.set("leftHeader", !0), t.g.vis.set("metaGaps", !0), t.g.vis.set("metaIdentity", !0), t.g.vis.set("metaLinks", !0)
                }), s.removeAllChilds(this.el), this.el.appendChild(this.buildDOM()), this
            },
            _addVisEl: function(t) {
                var e = this,
                    n = {};
                if (this.g.vis.get(t.id)) {
                    var r = "Hide ";
                    n.color = "red"
                } else r = "Show ", n.color = "green";
                return this.addNode(r + t.name, function() {
                    return e.g.vis.set(t.id, !e.g.vis.get(t.id))
                }, {
                    style: n
                })
            },
            getVisElements: function() {
                var t = [];
                return t.push({
                    name: "residues indices",
                    id: "markers"
                }), t.push({
                    name: "ID/Label",
                    id: "labels"
                }), t.push({
                    name: "meta info (Gaps/Ident)",
                    id: "metacell"
                }), t.push({
                    name: "overview panel",
                    id: "overviewbox"
                }), t.push({
                    name: "sequence logo",
                    id: "seqlogo"
                }), t.push({
                    name: "gap weights",
                    id: "gapHeader"
                }), t.push({
                    name: "conservation weights",
                    id: "conserv"
                }), t.push({
                    name: "scale slider",
                    id: "scaleslider"
                }), t.push({
                    name: "Label",
                    id: "labelName"
                }), t.push({
                    name: "ID",
                    id: "labelId"
                }), t.push({
                    name: "gaps %",
                    id: "metaGaps"
                }), t.push({
                    name: "identity score",
                    id: "metaIdentity"
                }), t
            }
        });
    e["default"] = u
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(35),
        o = r(i),
        s = n(26),
        u = r(s),
        a = n(27),
        l = r(a),
        c = n(28),
        f = r(c),
        h = n(29),
        d = r(h),
        p = n(9),
        g = r(p),
        v = n(30),
        m = r(v),
        y = n(32),
        _ = r(y),
        b = n(31),
        x = r(b),
        w = n(33),
        S = r(w),
        k = n(25),
        j = r(k),
        O = n(130),
        E = r(O),
        M = n(17),
        z = r(M),
        A = n(40),
        C = r(A),
        T = n(19),
        I = r(T),
        N = n(4),
        L = n(13),
        R = n(48),
        q = n(5),
        F = N.extend({
            initialize: function(t) {
                var e = this;
                if ("undefined" != typeof t && null !== t || (t = {}), null == t.colorscheme && (t.colorscheme = {}), null == t.columns && (t.columns = {}), null == t.conf && (t.conf = {}), null == t.vis && (t.vis = {}), null == t.visorder && (t.visorder = {}), null == t.zoomer && (t.zoomer = {}), null == t.conserv && (t.conserv = {}), null == t.scale && (t.scale = {}), this.g = L.mixin({}), this.seqs = this.g.seqs = new o["default"](t.seqs, this.g), this.g.config = new f["default"](t.conf), this.g["package"] = new d["default"](this.g), this.g.selcol = new g["default"]([], {
                        g: this.g
                    }), this.g.user = new m["default"], this.g.vis = new _["default"](t.vis, {
                        model: this.seqs
                    }), this.g.visorder = new x["default"](t.visorder), this.g.zoomer = new S["default"](t.zoomer, {
                        g: this.g,
                        model: this.seqs
                    }), this.g.scale = new j["default"](t.scale, {
                        g: this.g
                    }), this.g.conservationConfig = t.conserv, "localhost" === window.location.hostname && this.g.config.set("debug", !0), this._loadSeqs(t), this.u = {}, this.u.file = new z["default"](this), this.u.proxy = new I["default"]({
                        g: this.g
                    }), this.u.tree = new C["default"](this), this.g.config.get("eventBus") === !0 && this.startEventBus(), this.g.config.get("dropImport")) {
                    var n = {
                        dragover: this.dragOver,
                        drop: this.dropFile
                    };
                    this.delegateEvents(n)
                }
                return t.importURL && this.u.file.importURL(t.importURL, function() {
                    return e.render()
                }), t.bootstrapMenu && (t.menu && (this.menuConfig = t.menu), this.g.config.set("bootstrapMenu", !0)), this.draw(), this.m()
            },
            _loadSeqs: function(t) {
                var e = this.seqs.pluck("seq");
                return this.g.stats = new R(this.seqs, {
                    useGaps: !0
                }), this.g.stats.alphabetSize = this.g.config.get("alphabetSize"), this.g.columns = new l["default"](t.columns, this.g.stats), this.g.colorscheme = new u["default"](t.colorscheme, e, this.g.stats), this.g.zoomer.setEl(this.el, this.seqs)
            },
            importURL: function() {
                return this.u.file.importURL.apply(this.u.file, arguments)
            },
            m: function P() {
                var P = {};
                return P.model = n(51), P.selection = n(8), P.selcol = n(9), P.view = n(2), P.boneView = n(4), this.m = P
            },
            draw: function() {
                var t = this;
                if (this.removeViews(), this.addView("stage", new E["default"]({
                        model: this.seqs,
                        g: this.g
                    })), this.$el.addClass("biojs_msa_div"), this.g.config.get("bootstrapMenu")) {
                    var e = document.createElement("div"),
                        n = document.createElement("div");
                    this.el.parentNode ? (this.el.parentNode.replaceChild(n, this.el), n.appendChild(e), n.appendChild(this.el)) : (n.appendChild(e), n.appendChild(this.el));
                    var r = {
                        el: e,
                        msa: this
                    };
                    this.menuConfig && (r.menu = this.menuConfig), new msa.menu.defaultmenu(r).render()
                }
                return q(window).on("resize", function(e) {
                    var n = function() {
                        return this.g.zoomer.autoResize()
                    };
                    return setTimeout(n.bind(t), 5)
                })
            },
            dragOver: function(t) {
                return t.preventDefault(), t.target.className = "hover", !1
            },
            dropFile: function(t) {
                t.preventDefault();
                var e = t.target.files || t.dataTransfer.files;
                return this.u.file.importFiles(e), !1
            },
            startEventBus: function() {
                var t = this,
                    e = ["config", "columns", "colorscheme", "selcol", "vis", "visorder", "zoomer"];
                return function() {
                    for (var n, r = [], i = 0; i < e.length; i++) n = e[i], r.push(t._proxyToG(n));
                    return r
                }()
            },
            _proxyToG: function(t) {
                return this.listenTo(this.g[t], "all", function(e, n, r, i) {
                    if ("change" !== e) return "undefined" != typeof i && null !== i ? this.g.trigger(t + ":" + e, r, n, i) : this.g.trigger(t + ":" + e, r, n)
                })
            },
            render: function() {
                return void 0 === this.seqs || 0 === this.seqs.length, this.renderSubviews(), this.g.vis.set("loaded", !0), this
            }
        });
    e["default"] = F
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(134),
        o = r(i),
        s = n(143),
        u = r(s),
        a = n(4),
        l = a.extend({
            initialize: function(t) {
                this.g = t.g;
                var e = new u["default"]({
                    model: this.model,
                    g: this.g
                });
                if (e.ordering = -1, this.addView("labelblock", e), this.g.vis.get("sequences")) {
                    var n = new o["default"]({
                        model: this.model,
                        g: this.g
                    });
                    n.ordering = 0, this.addView("seqblock", n)
                }
                return this.listenTo(this.g.zoomer, "change:alignmentHeight", this.adjustHeight), this.listenTo(this.g.zoomer, "change:alignmentWidth", this.adjustWidth), this.listenTo(this.g.columns, "change:hidden", this.adjustHeight), this
            },
            render: function() {
                return this.renderSubviews(), this.el.className = "biojs_msa_albody", this.el.style.whiteSpace = "nowrap", this.adjustHeight(), this.adjustWidth(), this
            },
            adjustHeight: function() {
                return "auto" === this.g.zoomer.get("alignmentHeight") ? this.el.style.height = this.g.zoomer.get("rowHeight") * this.model.length + 5 : this.el.style.height = this.g.zoomer.get("alignmentHeight")
            },
            adjustWidth: function() {
                return this.el.style.width = this.getWidth()
            },
            getWidth: function() {
                var t = 0;
                return t += this.g.zoomer.getLeftBlockWidth(), this.g.vis.get("sequences") && (t += this.g.zoomer.get("alignmentWidth")), t
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(8),
        i = n(2),
        o = n(47),
        s = n(5),
        u = i.extend({
            className: "biojs_msa_overviewbox",
            tagName: "canvas",
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.zoomer, "change:boxRectWidth change:boxRectHeight change:overviewboxPaddingTop", this.rerender), this.listenTo(this.g.selcol, "add reset change", this.rerender), this.listenTo(this.g.columns, "change:hidden", this.rerender), this.listenTo(this.g.colorscheme, "change:showLowerCase", this.rerender), this.listenTo(this.model, "change", _.debounce(this.rerender, 5)), this.color = this.g.colorscheme.getSelectedScheme(), this.listenTo(this.g.colorscheme, "change:scheme", function() {
                    return this.color = this.g.colorscheme.getSelectedScheme(), this.rerender()
                }), this.dragStart = []
            },
            events: {
                click: "_onclick",
                mousedown: "_onmousedown"
            },
            rerender: function() {
                if (!this.g.config.get("manualRendering")) return this.render()
            },
            render: function() {
                this._createCanvas(), this.el.textContent = "overview", this.el.style.marginTop = this.g.zoomer.get("overviewboxPaddingTop"), this.ctx.fillStyle = "#999999", this.ctx.fillRect(0, 0, this.el.width, this.el.height);
                for (var t = this.g.zoomer.get("boxRectWidth"), e = this.g.zoomer.get("boxRectHeight"), n = this.g.columns.get("hidden"), r = this.g.colorscheme.get("showLowerCase"), i = -e, o = this.model.length, s = 0; s < o; s++)
                    if (this.model.at(s)) {
                        var u = this.model.at(s).get("seq"),
                            a = 0;
                        if (i += e, this.model.at(s).get("hidden")) this.ctx.fillStyle = "grey", this.ctx.fillRect(0, i, u.length * t, e);
                        else
                            for (var l = 0; l < u.length; l++) {
                                var c = u[l];
                                r && (c = c.toUpperCase());
                                var f = this.color.getColor(c, {
                                    pos: l
                                });
                                n.indexOf(l) >= 0 && (f = "grey"), "undefined" != typeof f && null !== f && (this.ctx.fillStyle = f, this.ctx.fillRect(a, i, t, e)), a += t
                            }
                    } return this._drawSelection()
            },
            _drawSelection: function() {
                var t = this;
                if (!(this.dragStart.length > 0) || this.prolongSelection) {
                    var e = this.g.zoomer.get("boxRectWidth"),
                        n = this.g.zoomer.get("boxRectHeight"),
                        r = n * this.model.length;
                    this.ctx.fillStyle = "#666666", this.ctx.globalAlpha = .9;
                    for (var i = this.g.selcol.length, o = function(i) {
                            var o = t.g.selcol.at(i);
                            if (!o) return "continue";
                            var s = void 0,
                                u = void 0;
                            "column" === o.get("type") ? t.ctx.fillRect(e * o.get("xStart"), 0, e * (o.get("xEnd") - o.get("xStart") + 1), r) : "row" === o.get("type") ? (s = t.model.filter(function(t) {
                                return t.get("id") === o.get("seqId")
                            })[0], u = t.model.indexOf(s), t.ctx.fillRect(0, n * u, e * s.get("seq").length, n)) : "pos" === o.get("type") && (s = t.model.filter(function(t) {
                                return t.get("id") === o.get("seqId")
                            })[0], u = t.model.indexOf(s), t.ctx.fillRect(e * o.get("xStart"), n * u, e * (o.get("xEnd") - o.get("xStart") + 1), n))
                        }, s = 0; s < i; s++) o(s);
                    return this.ctx.globalAlpha = 1
                }
            },
            _onclick: function(t) {
                return this.g.trigger("meta:click", {
                    seqId: this.model.get("id", {
                        evt: t
                    })
                })
            },
            _onmousemove: function(t) {
                if (0 !== this.dragStart.length) {
                    this.render(), this.ctx.fillStyle = "#666666", this.ctx.globalAlpha = .9;
                    var e = this._calcSelection(o.abs(t));
                    return this.ctx.fillRect(e[0][0], e[1][0], e[0][1] - e[0][0], e[1][1] - e[1][0]), t.preventDefault(), t.stopPropagation()
                }
            },
            _onmousedown: function(t) {
                var e = this;
                return this.dragStart = o.abs(t), this.dragStartRel = o.rel(t), t.ctrlKey || t.metaKey ? this.prolongSelection = !0 : this.prolongSelection = !1, s(document.body).on("mousemove.overmove", function(t) {
                    return e._onmousemove(t)
                }), s(document.body).on("mouseup.overup", function(t) {
                    return e._onmouseup(t)
                }), this.dragStart
            },
            _calcSelection: function(t) {
                for (var e = [t[0] - this.dragStart[0], t[1] - this.dragStart[1]], n = 0; n <= 1; n++) e[n] = this.dragStartRel[n] + e[n];
                for (var r = [
                        [this.dragStartRel[0], e[0]],
                        [this.dragStartRel[1], e[1]]
                    ], i = 0; i <= 1; i++) r[i][1] < r[i][0] && (r[i] = [r[i][1], r[i][0]]), r[i][0] = Math.max(r[i][0], 0);
                return r
            },
            _endSelection: function(t) {
                if (s(document.body).off(".overmove"), s(document.body).off(".overup"), 0 !== this.dragStart.length) {
                    for (var e = this._calcSelection(t), n = 0; n <= 1; n++) e[0][n] = Math.floor(e[0][n] / this.g.zoomer.get("boxRectWidth"));
                    for (var n = 0; n <= 1; n++) e[1][n] = Math.floor(e[1][n] / this.g.zoomer.get("boxRectHeight"));
                    e[0][1] = Math.min(this.model.getMaxLength() - 1, e[0][1]), e[1][1] = Math.min(this.model.length - 1, e[1][1]);
                    for (var i = [], o = e[1][0]; o <= e[1][1]; o++) {
                        var u = {
                            seqId: this.model.at(o).get("id"),
                            xStart: e[0][0],
                            xEnd: e[0][1]
                        };
                        i.push(new r.possel(u))
                    }
                    return this.dragStart = [], this.prolongSelection ? this.g.selcol.add(i) : this.g.selcol.reset(i), this.g.zoomer.setLeftOffset(e[0][0]), this.g.zoomer.setTopOffset(e[1][0])
                }
            },
            _onmouseup: function(t) {
                return this._endSelection(o.abs(t))
            },
            _onmouseout: function(t) {
                return this._endSelection(o.abs(t))
            },
            _createCanvas: function() {
                var t = this.g.zoomer.get("boxRectWidth"),
                    e = this.g.zoomer.get("boxRectHeight");
                return this.el.height = this.model.length * e, this.el.width = this.model.getMaxLength() * t, this.ctx = this.el.getContext("2d"), this.el.style.overflow = "auto", this.el.style.cursor = "crosshair"
            }
        });
    e["default"] = u
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(3),
        i = n(2),
        o = n(5),
        s = i.extend({
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.zoomer, "change:columnWidth", this.render), this.toggleClass = "msa-hide", this.isVisible = !0, this
            },
            attributes: {
                "class": "biojs_msa_scale"
            },
            events: {
                "change input": "updateSlider",
                "click button.msa-btn-close": "hide",
                "click button.msa-btn-open": "show",
                "click button[data-action]": "clickButton"
            },
            template: (0, r.template)('\t<div class="msa-scale-minimised">\t  <button class="btn msa-btn msa-btn-open">Zoom</button>\t</div>\t<div class="msa-scale-maximised">\t  <button class="btn msa-btn msa-btn-close" style="float:right">&times; close</button>\t  <div>\t  <input type="range" \t    data-provide="slider" \t    min="<%= min %>" \t    max="<%= max %>" \t    step="<%= step %>" \t    value="<%= value %>" \t  >\t  </div>\t  <div class="btngroup msa-btngroup">\t    <button class="btn msa-btn" data-action="smaller"><span class="glyphicon-zoom-out"></span>-</button>\t    <button class="btn msa-btn" data-action="bigger"><span class="glyphicon-zoom-in"></span>+</button>\t    <button class="btn msa-btn" data-action="reset"><span class="glyphicon-repeat"></span>reset</button>\t  </div>\t</div>\t'),
            render: function() {
                var t = this.model.getSizeRange(),
                    e = {
                        value: this.model.getSize(),
                        min: t[0],
                        max: t[1],
                        step: this.model.step || 1
                    };
                return this.$el.html(this.template(e)), this.isVisible ? this.show() : this.hide(), this
            },
            updateSlider: function(t) {
                var e = t.target,
                    n = parseInt(o(e).val());
                this.model.setSize(n)
            },
            clickButton: function(t) {
                var e = t.target,
                    n = o(e).data("action");
                return this.model[n], "function" == typeof this.model[n] && this.model[n](), this
            },
            hide: function() {
                this.isVisible = !1, this.$el.find(".msa-scale-minimised").removeClass(this.toggleClass), this.$el.find(".msa-scale-maximised").addClass(this.toggleClass)
            },
            show: function() {
                this.isVisible = !1, this.$el.find(".msa-scale-minimised").addClass(this.toggleClass), this.$el.find(".msa-scale-maximised").removeClass(this.toggleClass)
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(8),
        i = n(4),
        o = n(14),
        s = (n(7), i.extend({
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.user, "change:searchText", function(t, e) {
                    return this.search(e), this.render()
                }), this.sel = [], this.selPos = 0
            },
            events: {
                scroll: "_sendScrollEvent"
            },
            render: function() {
                this.renderSubviews(), this.el.className = "biojs_msa_searchresult";
                var t = this.g.user.get("searchText");
                return "undefined" != typeof t && null !== t && t.length > 0 && (0 === this.sel.length ? this.el.textContent = "no selection found" : (this.resultBox = o.mk("div"), this.resultBox.className = "biojs_msa_searchresult_ovbox", this.updateResult(), this.el.appendChild(this.resultBox), this.el.appendChild(this.buildBtns()))), this
            },
            updateResult: function() {
                var t = "search pattern: " + this.g.user.get("searchText");
                t += ", selection: " + (this.selPos + 1);
                var e = this.sel[this.selPos];
                return t += " (", t += e.get("xStart") + " - " + e.get("xEnd"), t += ", id: " + e.get("seqId"), t += ")", this.resultBox.textContent = t
            },
            buildBtns: function() {
                var t = this,
                    e = o.mk("button");
                e.textContent = "Prev", e.addEventListener("click", function() {
                    return t.moveSel(-1)
                });
                var n = o.mk("button");
                n.textContent = "Next", n.addEventListener("click", function() {
                    return t.moveSel(1)
                });
                var r = o.mk("button");
                r.textContent = "All", r.addEventListener("click", function() {
                    return t.g.selcol.reset(t.sel)
                });
                var i = o.mk("div");
                return i.appendChild(e), i.appendChild(n), i.appendChild(r), i.className = "biojs_msa_searchresult_row", i
            },
            moveSel: function(t) {
                var e = this.selPos + t;
                return e < 0 || e >= this.sel.length ? -1 : (this.focus(e), this.selPos = e, this.updateResult())
            },
            focus: function(t) {
                var e = this.sel[t],
                    n = e.get("xStart");
                return this.g.zoomer.setLeftOffset(n), this.g.selcol.reset([e])
            },
            search: function u(t) {
                var e, u = new RegExp(t, "gi"),
                    n = [],
                    i = e = 100042;
                return this.model.each(function(t) {
                    var e = t.get("seq");
                    return function() {
                        for (var o, s = []; o = u.exec(e);) {
                            var a = o.index,
                                l = {
                                    xStart: a,
                                    xEnd: a + o[0].length - 1,
                                    seqId: t.get("id")
                                };
                            n.push(new r.possel(l)), s.push(i = Math.min(a, i))
                        }
                        return s
                    }()
                }), this.g.selcol.reset(n), i === e && (i = 0), this.g.zoomer.setLeftOffset(i), this.sel = n
            }
        }));
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(3),
        o = n(126),
        s = r(o),
        u = n(138),
        a = r(u),
        l = n(127),
        c = r(l),
        f = n(129),
        h = r(f),
        d = n(128),
        p = r(d),
        g = n(4),
        v = g.extend({
            initialize: function(t) {
                return this.g = t.g, this.draw(), this.listenTo(this.g.stats, "reset", function() {
                    return this.rerender()
                }), this.listenTo(this.model, "change:hidden", (0, i.debounce)(this.rerender, 10)), this.listenTo(this.model, "sort", this.rerender), this.listenTo(this.model, "add", function() {
                    return console.log("seq add")
                }), this.listenTo(this.g.vis, "change:sequences", this.rerender), this.listenTo(this.g.vis, "change:overviewbox", this.rerender), this.listenTo(this.g.visorder, "change", this.rerender), this.listenTo(this.g.zoomer, "change:columnWidth", this.rerender), this.listenTo(this.g.vis, "change:scaleslider", this.rerender), this
            },
            draw: function() {
                if (this.removeViews(), this.g.vis.get("overviewbox")) {
                    var t = new c["default"]({
                        model: this.model,
                        g: this.g
                    });
                    t.ordering = this.g.visorder.get("overviewBox"), this.addView("overviewBox", t)
                }
                var e = new a["default"]({
                    model: this.model,
                    g: this.g
                });
                e.ordering = this.g.visorder.get("headerBox"), this.addView("headerBox", e);
                var n = new h["default"]({
                    model: this.model,
                    g: this.g
                });
                n.ordering = this.g.visorder.get("searchBox"), this.addView("searchbox", n);
                var r = new s["default"]({
                    model: this.model,
                    g: this.g
                });
                if (r.ordering = this.g.visorder.get("alignmentBody"), this.addView("body", r), this.g.vis.get("scaleslider")) {
                    var i = new p["default"]({
                        model: this.g.scale,
                        g: this.g
                    });
                    i.ordering = this.g.visorder.get("scaleSlider"), this.addView("scaleSlider", i)
                }
                return this
            },
            render: function(t) {
                return this.renderSubviews(), this.el.className = "biojs_msa_stage", this
            },
            rerender: function() {
                if (!this.g.config.get("manualRendering")) return this.draw(), this.render()
            }
        });
    e["default"] = v
}, function(t, e, n) {
    "use strict";

    function r(t, e) {
        if (!(t instanceof e)) throw new TypeError("Cannot call a class as a function")
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = function() {
            function t(t, e) {
                for (var n = 0; n < e.length; n++) {
                    var r = e[n];
                    r.enumerable = r.enumerable || !1, r.configurable = !0, "value" in r && (r.writable = !0), Object.defineProperty(t, r.key, r)
                }
            }
            return function(e, n, r) {
                return n && t(e.prototype, n), r && t(e, r), e
            }
        }(),
        o = (n(13), function() {
            function t(e) {
                r(this, t), this.g = e, this.cache = {}, this.cacheHeight = 0, this.cacheWidth = 0
            }
            return i(t, [{
                key: "getFontTile",
                value: function(t, e, n) {
                    return e === this.cacheWidth && n === this.cacheHeight || (this.cacheHeight = n, this.cacheWidth = e, this.cache = {}), void 0 === this.cache[t] && this.createTile(t, e, n), this.cache[t]
                }
            }, {
                key: "createTile",
                value: function(t, e, n) {
                    var r = this.cache[t] = document.createElement("canvas");
                    return r.width = e, r.height = n, this.ctx = r.getContext("2d"), this.ctx.font = this.g.zoomer.get("residueFont") + "px mono", this.ctx.textBaseline = "middle", this.ctx.textAlign = "center", this.ctx.fillText(t, e / 2, n / 2, e)
                }
            }]), t
        }());
    e["default"] = o
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(3),
        i = n(13),
        o = {
            setMaxScrollHeight: function() {
                return this.maxScrollHeight = this.g.zoomer.getMaxAlignmentHeight() - this.g.zoomer.get("alignmentHeight")
            },
            setMaxScrollWidth: function() {
                return this.maxScrollWidth = this.g.zoomer.getMaxAlignmentWidth() - this.g.zoomer.getAlignmentWidth()
            }
        },
        s = function(t, e) {
            return this.g = t, this.model = e, this.maxScrollWidth = 0, this.maxScrollHeight = 0, this.setMaxScrollHeight(), this.setMaxScrollWidth(), this.listenTo(this.g.zoomer, "change:rowHeight", this.setMaxScrollHeight), this.listenTo(this.g.zoomer, "change:columnWidth", this.setMaxScrollWidth), this.listenTo(this.g.zoomer, "change:alignmentWidth", this.setMaxScrollWidth), this.listenTo(this.g.zoomer, "change:alignmentHeight", this.setMaxScrollHeight), this.listenTo(this.model, "add change reset", function() {
                return this.setMaxScrollHeight(), this.setMaxScrollWidth()
            }, this), this
        };
    (0, r.extend)(s.prototype, o), i.mixin(s.prototype), e["default"] = s
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(3),
        i = function(t, e) {
            return this.g = t, this.ctx = e, this
        };
    _.extend(i.prototype, {
        _getSelection: function(t) {
            var e = t.get("seq").length,
                n = [],
                i = this.g.selcol.getSelForRow(t.get("id")),
                o = (0, r.find)(i, function(t) {
                    return "row" === t.get("type")
                });
            if ("undefined" != typeof o && null !== o)
                for (var s = e - 1, u = 0; u <= s; u++) n.push(u);
            else if (i.length > 0)
                for (var a, l = 0; l < i.length; l++) {
                    a = i[l];
                    for (var c = a.get("xStart"), f = a.get("xEnd"), h = c; h <= f; h++) n.push(h)
                }
            return n
        },
        _appendSelection: function(t) {
            var e = this,
                n = t.model.get("seq"),
                r = this._getSelection(t.model),
                i = this._getPrevNextSelection(t.model),
                o = i[0],
                s = i[1];
            if (this.g.zoomer.get("columnWidth"), this.g.zoomer.get("rowHeight"), 0 !== r.length) {
                var u = 0;
                return function() {
                    for (var i = [], a = n.length - 1, l = function(n) {
                            i.push(function() {
                                if (t.hidden.indexOf(n) >= 0) return u++;
                                var i = n - u;
                                return r.indexOf(n) >= 0 && (0 === i || r.indexOf(n - 1) < 0) ? e._renderSelection({
                                    n: n,
                                    k: i,
                                    selection: r,
                                    mPrevSel: o,
                                    mNextSel: s,
                                    xZero: t.xZero,
                                    yZero: t.yZero,
                                    model: t.model
                                }) : void 0
                            }())
                        }, c = 0; 0 < a ? c <= a : c >= a; 0 < a ? c++ : c--) l(c);
                    return i
                }()
            }
        },
        _renderSelection: function(t) {
            for (var e = t.xZero, n = t.yZero, r = t.n, i = t.k, o = t.selection, s = t.mPrevSel, u = t.mNextSel, a = 0, l = t.model.get("seq").length - 1, c = r;
                (r < l ? c <= l : c >= l) && o.indexOf(c) >= 0; r < l ? c++ : c--) a++;
            var f = this.g.zoomer.get("columnWidth"),
                h = this.g.zoomer.get("rowHeight"),
                d = f * a + 1,
                p = this.g.columns.get("hidden");
            this.ctx.beginPath();
            var g = this.ctx.lineWidth;
            this.ctx.lineWidth = 3;
            var v = this.ctx.strokeStyle;
            this.ctx.strokeStyle = "#FF0000", e += i * f;
            for (var m = 0, y = a - 1, _ = 0; 0 < y ? _ <= y : _ >= y; 0 < y ? _++ : _--) {
                var b = r + _;
                p.indexOf(b) >= 0 || ("undefined" != typeof s && null !== s && s.indexOf(b) >= 0 || (this.ctx.moveTo(e + m, n),
                    this.ctx.lineTo(m + f + e, n)), "undefined" != typeof u && null !== u && u.indexOf(b) >= 0 || (this.ctx.moveTo(m + e, h + n), this.ctx.lineTo(m + f + e, h + n)), m += f)
            }
            return this.ctx.moveTo(e, n), this.ctx.lineTo(e, h + n), this.ctx.moveTo(e + d, n), this.ctx.lineTo(e + d, h + n), this.ctx.stroke(), this.ctx.strokeStyle = v, this.ctx.lineWidth = g
        },
        _getPrevNextSelection: function(t) {
            var e = t.collection.prev(t),
                n = t.collection.next(t),
                r = void 0,
                i = void 0;
            return "undefined" != typeof e && null !== e && (r = this._getSelection(e)), "undefined" != typeof n && null !== n && (i = this._getSelection(n)), [r, i]
        }
    }), e["default"] = i
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = function() {
            function t(t, e) {
                var n = [],
                    r = !0,
                    i = !1,
                    o = void 0;
                try {
                    for (var s, u = t[Symbol.iterator](); !(r = (s = u.next()).done) && (n.push(s.value), !e || n.length !== e); r = !0);
                } catch (a) {
                    i = !0, o = a
                } finally {
                    try {
                        !r && u["return"] && u["return"]()
                    } finally {
                        if (i) throw o
                    }
                }
                return n
            }
            return function(e, n) {
                if (Array.isArray(e)) return e;
                if (Symbol.iterator in Object(e)) return t(e, n);
                throw new TypeError("Invalid attempt to destructure non-iterable instance")
            }
        }(),
        o = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        },
        s = n(3),
        u = n(131),
        a = r(u),
        l = n(133),
        c = r(l),
        f = n(135),
        h = r(f),
        d = n(132),
        p = r(d),
        g = n(4),
        v = n(47),
        m = n(5),
        y = g.extend({
            tagName: "canvas",
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.zoomer, "change:_alignmentScrollLeft change:_alignmentScrollTop", function(t, e, n) {
                    if (null == ("undefined" != typeof n && null !== n ? n.origin : void 0) || "canvasseq" !== n.origin) return this.render()
                }), this.listenTo(this.g.columns, "change:hidden", this.render), this.listenTo(this.g.zoomer, "change:alignmentWidth change:alignmentHeight", this.render), this.listenTo(this.g.colorscheme, "change", this.render), this.listenTo(this.g.selcol, "reset add", this.render), this.listenTo(this.model, "reset add", this.render), this.el.style.display = "inline-block", this.el.style.overflowX = "hidden", this.el.style.overflowY = "hidden", this.el.className = "biojs_msa_seqblock", this.ctx = this.el.getContext("2d"), this.cache = new a["default"](this.g), this.coordsCache = new p["default"](this.g, this.model), this.listenTo(this.g.zoomer, "change:residueFont", function() {
                    return this.cache = new a["default"](this.g), this.render()
                }), this.sel = new c["default"](this.g, this.ctx), this._setColor(), this.throttleTime = 0, this.throttleCounts = 0, null != document.documentElement.style.webkitAppearance ? this.throttledDraw = function() {
                    var t = +new Date;
                    if (this.draw(), this.throttleTime += +new Date - t, this.throttleCounts++, this.throttleCounts > 15) return Math.ceil(this.throttleTime / this.throttleCounts), this.throttledDraw = this.draw
                } : this.throttledDraw = (0, s.throttle)(this.throttledDraw, 30), this.manageEvents()
            },
            throttledDraw: function() {
                var t = +new Date;
                if (this.draw(), this.throttleTime += +new Date - t, this.throttleCounts++, this.throttleCounts > 15) {
                    var e = Math.ceil(this.throttleTime / this.throttleCounts);
                    return e *= 1.2, e = Math.max(20, e), this.throttledDraw = _.throttle(this.draw, e)
                }
            },
            manageEvents: function() {
                var t = {};
                return t.mousedown = "_onmousedown", t.touchstart = "_ontouchstart", this.g.config.get("registerMouseClicks") && (t.dblclick = "_onclick"), this.g.config.get("registerMouseHover") && (t.mousein = "_onmousein", t.mouseout = "_onmouseout"), t.mousewheel = "_onmousewheel", t.DOMMouseScroll = "_onmousewheel", this.delegateEvents(t), this.listenTo(this.g.config, "change:registerMouseHover", this.manageEvents), this.listenTo(this.g.config, "change:registerMouseClick", this.manageEvents), this.dragStart = []
            },
            _setColor: function() {
                return this.color = this.g.colorscheme.getSelectedScheme()
            },
            draw: function() {
                if (this.el.width = this.el.width, null != this.seqDrawer && this.model.length > 0) return this.seqDrawer.drawLetters(), this.seqDrawer.drawRows(this.sel._appendSelection, this.sel), this.seqDrawer.drawRows(this.drawFeatures, this)
            },
            drawFeatures: function(t) {
                var e = this,
                    n = this.g.zoomer.get("columnWidth"),
                    r = this.g.zoomer.get("rowHeight");
                if (t.model.attributes.height > 1) {
                    var i = function() {
                        var i = e.ctx;
                        return t.model.attributes.features.each(function(e) {
                            i.fillStyle = e.attributes.fillColor || "red";
                            var o = e.attributes.xEnd - e.attributes.xStart + 1,
                                s = (e.attributes.row + 1) * r;
                            return i.fillRect(e.attributes.xStart * n + t.xZero, s + t.yZero, n * o, r)
                        }), i.fillStyle = "black", i.font = e.g.zoomer.get("residueFont") + "px mono", i.textBaseline = "middle", i.textAlign = "center", {
                            v: t.model.attributes.features.each(function(e) {
                                var o = e.attributes.xEnd - e.attributes.xStart + 1,
                                    s = (e.attributes.row + 1) * r;
                                return i.fillText(e.attributes.text, t.xZero + e.attributes.xStart * n + o / 2 * n, t.yZero + .5 * r + s)
                            })
                        }
                    }();
                    if ("object" === ("undefined" == typeof i ? "undefined" : o(i))) return i.v
                }
            },
            render: function() {
                return this.el.setAttribute("height", this.g.zoomer.get("alignmentHeight") + "px"), this.el.setAttribute("width", this.g.zoomer.getAlignmentWidth() + "px"), this.g.zoomer._checkScrolling(this._checkScrolling([this.g.zoomer.get("_alignmentScrollLeft"), this.g.zoomer.get("_alignmentScrollTop")]), {
                    header: "canvasseq"
                }), this._setColor(), this.seqDrawer = new h["default"](this.g, this.ctx, this.model, {
                    width: this.el.width,
                    height: this.el.height,
                    color: this.color,
                    cache: this.cache
                }), this.throttledDraw(), this
            },
            _onmousemove: function(t, e) {
                if (0 !== this.dragStart.length) {
                    var n = v.abs(t),
                        r = [n[0] - this.dragStart[0], n[1] - this.dragStart[1]],
                        i = this.g.zoomer.get("canvasEventScale");
                    e && (i = 3);
                    for (var o = 0; o <= 1; o++) r[o] = r[o] * i;
                    for (var s = [this.dragStartScroll[0] - r[0], this.dragStartScroll[1] - r[1]], u = 0; u <= 1; u++) s[u] = Math.round(s[u]);
                    var a = this._checkScrolling(s);
                    this.g.zoomer._checkScrolling(a, {
                        origin: "canvasseq"
                    });
                    for (var l = 0; l <= 1; l++) a[l] !== s[l] && (0 === a[l] ? (this.dragStart[l] = n[l], this.dragStartScroll[l] = 0) : this.dragStart[l] = n[l] - a[l]);
                    return this.throttledDraw(), null != t.preventDefault ? (t.preventDefault(), t.stopPropagation()) : void 0
                }
            },
            _ontouchmove: function(t) {
                return this._onmousemove(t.changedTouches[0], !0), t.preventDefault(), t.stopPropagation()
            },
            _onmousedown: function(t) {
                var e = this;
                return this.dragStart = v.abs(t), this.dragStartScroll = [this.g.zoomer.get("_alignmentScrollLeft"), this.g.zoomer.get("_alignmentScrollTop")], m(document.body).on("mousemove.overmove", function(t) {
                    return e._onmousemove(t)
                }), m(document.body).on("mouseup.overup", function() {
                    return e._cleanup()
                }), t.preventDefault()
            },
            _ontouchstart: function(t) {
                var e = this;
                return this.dragStart = v.abs(t.changedTouches[0]), this.dragStartScroll = [this.g.zoomer.get("_alignmentScrollLeft"), this.g.zoomer.get("_alignmentScrollTop")], m(document.body).on("touchmove.overtmove", function(t) {
                    return e._ontouchmove(t)
                }), m(document.body).on("touchend.overtend touchleave.overtleave touchcancel.overtcanel", function(t) {
                    return e._touchCleanup(t)
                })
            },
            _onmousewinout: function(t) {
                if (t.toElement === document.body.parentNode) return this._cleanup()
            },
            _cleanup: function() {
                return this.dragStart = [], m(document.body).off(".overmove"), m(document.body).off(".overup"), m(document.body).off(".overout")
            },
            _touchCleanup: function(t) {
                return t.changedTouches.length > 0 && this._onmousemove(t.changedTouches[0], !0), this.dragStart = [], m(document.body).off(".overtmove"), m(document.body).off(".overtend"), m(document.body).off(".overtleave"), m(document.body).off(".overtcancel")
            },
            _onmousewheel: function(t) {
                var e = v.wheelDelta(t);
                return this.g.zoomer.set("_alignmentScrollLeft", this.g.zoomer.get("_alignmentScrollLeft") + e[0]), this.g.zoomer.set("_alignmentScrollTop", this.g.zoomer.get("_alignmentScrollTop") + e[1]), t.preventDefault()
            },
            _onclick: function(t) {
                var e = this._getClickPos(t);
                return "undefined" != typeof e && null !== e && (null != e.feature ? this.g.trigger("feature:click", e) : this.g.trigger("residue:click", e)), this.throttledDraw()
            },
            _onmousein: function(t) {
                var e = this._getClickPos(t);
                return "undefined" != typeof e && null !== e && (null != e.feature ? this.g.trigger("feature:mousein", e) : this.g.trigger("residue:mousein", e)), this.throttledDraw()
            },
            _onmouseout: function(t) {
                var e = this._getClickPos(t);
                return "undefined" != typeof e && null !== e && (null != e.feature ? this.g.trigger("feature:mouseout", e) : this.g.trigger("residue:mouseout", e)), this.throttledDraw()
            },
            _getClickPos: function(t) {
                var e = v.rel(t);
                e[0] += this.g.zoomer.get("_alignmentScrollLeft");
                var n = Math.floor(e[0] / this.g.zoomer.get("columnWidth")),
                    r = this.seqDrawer._getSeqForYClick(e[1]),
                    o = i(r, 2),
                    s = o[0],
                    u = o[1];
                n += this.g.columns.calcHiddenColumns(n), s += this.model.calcHiddenSeqs(s), n = Math.max(0, n), s = Math.max(0, s);
                var a = this.model.at(s).get("id");
                if (!(u > 0)) return {
                    seqId: a,
                    rowPos: n,
                    evt: t
                };
                var l = this.model.at(s).get("features").getFeatureOnRow(u - 1, n);
                return 0 !== l.length ? {
                    seqId: a,
                    feature: l[0],
                    rowPos: n,
                    evt: t
                } : void 0
            },
            _checkScrolling: function(t) {
                for (var e = [this.coordsCache.maxScrollWidth, this.coordsCache.maxScrollHeight], n = 0; n <= 1; n++) t[n] > e[n] && (t[n] = e[n]), t[n] < 0 && (t[n] = 0);
                return t
            }
        });
    e["default"] = y
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = function() {
            function t(t, e) {
                var n = [],
                    r = !0,
                    i = !1,
                    o = void 0;
                try {
                    for (var s, u = t[Symbol.iterator](); !(r = (s = u.next()).done) && (n.push(s.value), !e || n.length !== e); r = !0);
                } catch (a) {
                    i = !0, o = a
                } finally {
                    try {
                        !r && u["return"] && u["return"]()
                    } finally {
                        if (i) throw o
                    }
                }
                return n
            }
            return function(e, n) {
                if (Array.isArray(e)) return e;
                if (Symbol.iterator in Object(e)) return t(e, n);
                throw new TypeError("Invalid attempt to destructure non-iterable instance")
            }
        }(),
        i = n(3),
        o = {
            updateConfig: function() {
                this.rectWidth = this.g.zoomer.get("columnWidth"), this.rectHeight = this.g.zoomer.get("rowHeight")
            },
            drawLetters: function() {
                return this.updateConfig(), this.ctx.globalAlpha = this.g.colorscheme.get("opacity"), this.drawSeqs(function(t) {
                    return this.drawSeq(t, this._drawRect)
                }), this.ctx.globalAlpha = 1, this.rectWidth >= this.g.zoomer.get("minLetterDrawSize") && this.drawSeqs(function(t) {
                    return this.drawSeq(t, this._drawLetter)
                }), this
            },
            drawSeqs: function(t, e) {
                var n = this.g.columns.get("hidden");
                e = e || this;
                for (var i = this.getStartSeq(), o = r(i, 2), s = o[0], u = o[1], a = s; a < this.model.length; a++) {
                    var l = this.model.at(a);
                    if (!l.get("hidden") && (t.call(e, {
                            model: l,
                            yPos: u,
                            y: a,
                            hidden: n
                        }), u += (l.attributes.height || 1) * this.rectHeight, u > this.height)) break
                }
            },
            drawRows: function(t, e) {
                return this.drawSeqs(function(n) {
                    return this.drawRow(n, t, e)
                })
            },
            drawRow: function(t, e, n) {
                var r = this.g.zoomer.get("columnWidth"),
                    i = Math.max(0, Math.abs(Math.ceil(-this.g.zoomer.get("_alignmentScrollLeft") / r))),
                    o = -Math.abs(-this.g.zoomer.get("_alignmentScrollLeft") % r),
                    s = o - i * r,
                    u = t.yPos;
                return e.call(n, {
                    model: t.model,
                    xZero: s,
                    yZero: u,
                    hidden: t.hidden
                })
            },
            getStartSeq: function() {
                for (var t = Math.max(0, Math.floor(this.g.zoomer.get("_alignmentScrollTop") / this.rectHeight)) + 1, e = 0, n = 0; e < t && n < this.model.length;) e += this.model.at(n).attributes.height || 1, n++;
                return [n - 1, -Math.max(0, this.g.zoomer.get("_alignmentScrollTop") - e * this.rectHeight + (this.model.at(n - 1).attributes.height || 1) * this.rectHeight)]
            },
            _getSeqForYClick: function(t) {
                for (var e = this.getStartSeq(), n = r(e, 2), i = n[0], o = n[1], s = o % this.rectHeight, u = Math.max(0, Math.floor((t - s) / this.rectHeight)) + 1, a = 0, l = i; a < u && l < this.model.length;) a += this.model.at(l).attributes.height || 1, l++;
                return [l - 1, Math.max(0, Math.floor(t / this.rectHeight) - a + (this.model.at(l - 1).get("height") || 1))]
            },
            drawSeq: function(t, e) {
                for (var n = t.model.get("seq"), r = t.yPos, i = this.rectWidth, o = this.rectHeight, s = Math.max(0, Math.abs(Math.ceil(-this.g.zoomer.get("_alignmentScrollLeft") / i))), u = -Math.abs(-this.g.zoomer.get("_alignmentScrollLeft") % i), a = {
                        rectWidth: i,
                        rectHeight: o,
                        yPos: r,
                        y: t.y
                    }, l = this.width, c = s; c < n.length; c++) {
                    var f = n[c];
                    if (f = f.toUpperCase(), a.x = c, a.c = f, a.xPos = u, t.hidden.indexOf(c) < 0 && (e(this, a), u += i, u > l)) break
                }
            },
            _drawRect: function(t, e) {
                var n = t.color.getColor(e.c, {
                    pos: e.x,
                    y: e.y
                });
                if ("undefined" != typeof n && null !== n) return t.ctx.fillStyle = n, t.ctx.fillRect(e.xPos, e.yPos, e.rectWidth, e.rectHeight)
            },
            _drawLetter: function(t, e) {
                return t.ctx.drawImage(t.cache.getFontTile(e.c, e.rectWidth, e.rectHeight), e.xPos, e.yPos, e.rectWidth, e.rectHeight)
            }
        },
        s = function(t, e, n, r) {
            return this.g = t, this.ctx = e, this.model = n, this.width = r.width, this.height = r.height, this.color = r.color, this.cache = r.cache, this.rectHeight = this.g.zoomer.get("rowHeight"), this.rectWidth = this.g.zoomer.get("columnWidth"), this
        };
    (0, i.extend)(s.prototype, o), e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        if (t && t.__esModule) return t;
        var e = {};
        if (null != t)
            for (var n in t) Object.prototype.hasOwnProperty.call(t, n) && (e[n] = t[n]);
        return e["default"] = t, e
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
            return typeof t
        } : function(t) {
            return t && "function" == typeof Symbol && t.constructor === Symbol ? "symbol" : typeof t
        },
        o = n(12),
        s = r(o),
        u = n(2),
        a = n(7),
        l = u.extend({
            className: "biojs_msa_conserv",
            initialize: function(t) {
                this.g = t.g, this.listenTo(this.g.zoomer, "change:stepSize change:labelWidth change:columnWidth", this.render), this.listenTo(this.g.vis, "change:labels change:metacell", this.render), this.listenTo(this.g.columns, "change:scaling", this.render), this.listenTo(this.g.stats, "reset", this.render);
                var e = _.extend({}, {
                    fillColor: ["#660", "#ff0"],
                    strokeColor: "#330",
                    maxHeight: 20,
                    rectStyler: function(t, e) {
                        return t
                    }
                }, this.g.conservationConfig);
                return this.fillColor = e.fillColor, this.strokeColor = e.strokeColor, this.maxHeight = e.maxHeight, this.rectStyler = e.rectStyler, this.manageEvents()
            },
            colorer: function c(t) {
                var e = this,
                    c = function() {
                        return "none"
                    };
                if ("string" == typeof t) c = function() {
                    return t
                };
                else if (Array.isArray(t)) {
                    2 != t.length && console.error("ERROR: colorRange array should have exactly two elements", t);
                    var n = "undefined" != typeof d3 && !!d3.scale,
                        r = "undefined" != typeof d3_scale;
                    n || r ? ! function() {
                        var r = n ? d3.scale.linear() : d3_scale.scaleLinear(),
                            i = r.domain([0, e.maxHeight]).range(t);
                        c = function(t) {
                            return i(t.height)
                        }
                    }() : (console.warn("providing a [min/max] range as input requires d3 to be included - only using the first color"), c = function(e) {
                        return t[0]
                    })
                } else console.warn("expected colorRange to be '#rgb' or ['#rgb', '#rgb']", t, "(" + ("undefined" == typeof t ? "undefined" : i(t)) + ")");
                return c
            },
            render: function() {
                var t = this.g.stats.scale(this.g.stats.conservation());
                a.removeAllChilds(this.el);
                var e = this.model.getMaxLength(),
                    n = this.g.zoomer.get("columnWidth"),
                    r = this.maxHeight,
                    i = n * (e - this.g.columns.get("hidden").length),
                    o = s.base({
                        height: r,
                        width: i
                    });
                o.style.display = "inline-block", o.style.cursor = "pointer";
                for (var u = (this.rectData, this.colorer(this.fillColor)), l = this.colorer(this.strokeColor), c = this.rectStyler, f = this.g.zoomer.get("stepSize"), h = this.g.columns.get("hidden"), d = 0, p = 0; p < e;)
                    if (h.indexOf(p) >= 0) p += f;
                    else {
                        i = n * f;
                        for (var g = 0, v = f - 1, m = 0; 0 < v ? m <= v : m >= v; 0 < v ? m++ : m--) g += t[p];
                        var y = r * (g / f),
                            _ = {
                                x: d,
                                y: r - y,
                                maxheight: r,
                                width: i - n / 4,
                                height: y,
                                rowPos: p
                            },
                            b = s.rect(_);
                        b.style.stroke = l(_), b.style.fill = u(_), "function" == typeof c && c(b, _), b.rowPos = p, o.appendChild(b), d += i, p += f
                    } return this.el.appendChild(o), this
            },
            _onclick: function(t) {
                var e = this,
                    n = t.target.rowPos,
                    r = this.g.zoomer.get("stepSize");
                return function() {
                    for (var i = [], o = r - 1, s = 0; 0 < o ? s <= o : s >= o; 0 < o ? s++ : s--) i.push(e.g.trigger("bar:click", {
                        rowPos: n + s,
                        evt: t
                    }));
                    return i
                }()
            },
            manageEvents: function() {
                var t = {};
                return this.g.config.get("registerMouseClicks") && (t.click = "_onclick"), this.g.config.get("registerMouseHover") && (t.mousein = "_onmousein", t.mouseout = "_onmouseout"), this.delegateEvents(t), this.listenTo(this.g.config, "change:registerMouseHover", this.manageEvents), this.listenTo(this.g.config, "change:registerMouseClick", this.manageEvents)
            },
            _onmousein: function(t) {
                var e = this.g.zoomer.get("stepSize" * t.rowPos);
                return this.g.trigger("bar:mousein", {
                    rowPos: e,
                    evt: t
                })
            },
            _onmouseout: function(t) {
                var e = this.g.zoomer.get("stepSize" * t.rowPos);
                return this.g.trigger("bar:mouseout", {
                    rowPos: e,
                    evt: t
                })
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";

    function r(t) {
        if (t && t.__esModule) return t;
        var e = {};
        if (null != t)
            for (var n in t) Object.prototype.hasOwnProperty.call(t, n) && (e[n] = t[n]);
        return e["default"] = t, e
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(12),
        o = r(i),
        s = n(2),
        u = n(7),
        a = s.extend({
            className: "biojs_msa_gapview",
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.zoomer, "change:stepSize change:labelWidth change:columnWidth", this.render), this.listenTo(this.g.vis, "change:labels change:metacell", this.render), this.listenTo(this.g.columns, "change:scaling", this.render), this.listenTo(this.model, "reset", this.render), this.manageEvents()
            },
            render: function() {
                var t = this.g.stats.gaps();
                u.removeAllChilds(this.el);
                var e = this.model.getMaxLength(),
                    n = this.g.zoomer.get("columnWidth"),
                    r = 20,
                    i = n * (e - this.g.columns.get("hidden").length),
                    s = o.base({
                        height: r,
                        width: i
                    });
                s.style.display = "inline-block", s.style.cursor = "pointer";
                for (var a = this.g.zoomer.get("stepSize"), l = this.g.columns.get("hidden"), c = 0, f = 0; f < e;)
                    if (l.indexOf(f) >= 0) f += a;
                    else {
                        i = n * a;
                        for (var h = 0, d = a - 1, p = 0; 0 < d ? p <= d : p >= d; 0 < d ? p++ : p--) h += t[f];
                        var g = r * (h / a),
                            v = o.rect({
                                x: c,
                                y: r - g,
                                width: i - n / 4,
                                height: g,
                                style: "stroke:red;stroke-width:1;"
                            });
                        v.rowPos = f, s.appendChild(v), c += i, f += a
                    } return this.el.appendChild(s), this
            },
            _onclick: function(t) {
                var e = this,
                    n = t.target.rowPos,
                    r = this.g.zoomer.get("stepSize");
                return function() {
                    for (var i = [], o = r - 1, s = 0; 0 < o ? s <= o : s >= o; 0 < o ? s++ : s--) i.push(e.g.trigger("gap:click", {
                        rowPos: n + s,
                        evt: t
                    }));
                    return i
                }()
            },
            manageEvents: function() {
                var t = {};
                return this.g.config.get("registerMouseClicks") && (t.click = "_onclick"), this.g.config.get("registerMouseHover") && (t.mousein = "_onmousein", t.mouseout = "_onmouseout"), this.delegateEvents(t), this.listenTo(this.g.config, "change:registerMouseHover", this.manageEvents), this.listenTo(this.g.config, "change:registerMouseClick", this.manageEvents)
            },
            _onmousein: function(t) {
                var e = this.g.zoomer.get("stepSize" * t.rowPos);
                return this.g.trigger("gap:mousein", {
                    rowPos: e,
                    evt: t
                })
            },
            _onmouseout: function(t) {
                var e = this.g.zoomer.get("stepSize" * t.rowPos);
                return this.g.trigger("gap:mouseout", {
                    rowPos: e,
                    evt: t
                })
            }
        });
    e["default"] = a
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(139),
        o = r(i),
        s = n(141),
        u = r(s),
        a = n(4),
        l = a.extend({
            initialize: function(t) {
                var e = this;
                return this.g = t.g, this.draw(), this.listenTo(this.g.vis, "change:labels change:metacell change:leftHeader", function() {
                    return e.draw(), e.render()
                })
            },
            draw: function() {
                if (this.removeViews(), this.g.vis.get("leftHeader") && (this.g.vis.get("labels") || this.g.vis.get("metacell"))) {
                    var t = new o["default"]({
                        model: this.model,
                        g: this.g
                    });
                    t.ordering = -50, this.addView("lHeader", t)
                }
                var e = new u["default"]({
                    model: this.model,
                    g: this.g
                });
                return e.ordering = 0, this.addView("rHeader", e)
            },
            render: function() {
                return this.renderSubviews(), this.el.className = "biojs_msa_header"
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(14),
        i = n(2),
        o = n(7),
        s = i.extend({
            className: "biojs_msa_headers",
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.vis, "change:metacell change:labels", this.render), this.listenTo(this.g.zoomer, "change:labelWidth change:metaWidth", this.render)
            },
            render: function() {
                o.removeAllChilds(this.el);
                var t = 0;
                return t += this.g.zoomer.getLeftBlockWidth(), this.el.style.width = t + "px", this.g.vis.get("labels") && this.el.appendChild(this.labelDOM()), this.g.vis.get("metacell") && this.el.appendChild(this.metaDOM()), this.el.style.display = "inline-block", this.el.style.fontSize = this.g.zoomer.get("markerFontsize"), this
            },
            labelDOM: function() {
                var t = r.mk("div");
                if (t.style.width = this.g.zoomer.getLabelWidth(), t.style.display = "inline-block", this.g.vis.get("labelCheckbox") && t.appendChild(this.addEl(".", 10)), this.g.vis.get("labelId") && t.appendChild(this.addEl("ID", this.g.zoomer.get("labelIdLength"))), this.g.vis.get("labelPartition") && t.appendChild(this.addEl("part", 15)), this.g.vis.get("labelName")) {
                    var e = this.addEl("Label");
                    t.appendChild(e)
                }
                return t
            },
            addEl: function(t, e) {
                var n = document.createElement("span");
                return n.textContent = t, "undefined" != typeof e && null !== e && (n.style.width = e + "px"), n.style.display = "inline-block", n
            },
            metaDOM: function() {
                var t = r.mk("div");
                return t.style.width = this.g.zoomer.getMetaWidth(), t.style.display = "inline-block", this.g.vis.get("metaGaps") && t.appendChild(this.addEl("Gaps", this.g.zoomer.get("metaGapWidth"))), this.g.vis.get("metaIdentity") && t.appendChild(this.addEl("Ident", this.g.zoomer.get("metaIdentWidth"))), t
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        if (t && t.__esModule) return t;
        var e = {};
        if (null != t)
            for (var n in t) Object.prototype.hasOwnProperty.call(t, n) && (e[n] = t[n]);
        return e["default"] = t, e
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(12),
        o = r(i),
        s = n(2),
        u = n(7),
        a = n(5),
        l = s.extend({
            className: "biojs_msa_marker",
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.zoomer, "change:stepSize change:labelWidth change:columnWidth change:markerStepSize change:markerFontsize", this.render), this.listenTo(this.g.vis, "change:labels change:metacell", this.render), this.manageEvents()
            },
            render: function() {
                u.removeAllChilds(this.el);
                var t = this.g.zoomer.get("markerFontsize"),
                    e = this.g.zoomer.get("columnWidth"),
                    n = this.g.zoomer.get("stepSize"),
                    r = this.g.zoomer.get("markerStepSize"),
                    i = this.g.columns.get("hidden");
                this.el.style.fontSize = t;
                for (var o = document.createElement("span"), s = 0, a = this.model.getMaxLength(), s = 0; s < a; s++)
                    if (i.indexOf(s) >= 0) this.markerHidden(l, s, n), s += n;
                    else {
                        var l = document.createElement("span");
                        l.className = "msa-col-header", l.style.width = e + "px", l.style.display = "inline-block", (s + 1) % r === 0 ? l.textContent = s + 1 : (s + 1) % n === 0 ? l.textContent = "." : l.textContent = " ", l.rowPos = s, o.appendChild(l)
                    } return this.el.appendChild(o), this
            },
            markerHidden: function(t, e, n) {
                for (var r = this, i = this.g.columns.get("hidden").slice(0), s = Math.max(0, e - n), u = !0, l = s; s < e ? l <= e : l >= e; s < e ? l++ : l--) u &= i.indexOf(l) >= 0;
                if (!u) {
                    for (var c = this.model.getMaxLength(), f = 0, h = -1, e = e;
                        (e < c ? e <= c : e >= c) && (h >= 0 || (h = i.indexOf(e)), i.indexOf(e) >= 0); e < c ? e++ : e--) f++;
                    var d = o.base({
                        height: 10,
                        width: 10
                    });
                    d.style.position = "relative";
                    var p = o.polygon({
                        points: "0,0 5,5 10,0",
                        style: "fill:lime;stroke:purple;stroke-width:1"
                    });
                    return a(p).on("click", function(t) {
                        return i.splice(h, f), r.g.columns.set("hidden", i)
                    }), d.appendChild(p), t.appendChild(d), d
                }
            },
            manageEvents: function() {
                var t = {};
                return this.g.config.get("registerMouseClicks") && (t.click = "_onclick"), this.g.config.get("registerMouseHover") && (t.mousein = "_onmousein", t.mouseout = "_onmouseout"), this.delegateEvents(t), this.listenTo(this.g.config, "change:registerMouseHover", this.manageEvents), this.listenTo(this.g.config, "change:registerMouseClick", this.manageEvents)
            },
            _onclick: function(t) {
                var e = t.target.rowPos,
                    n = this.g.zoomer.get("stepSize");
                return this.g.trigger("column:click", {
                    rowPos: e,
                    stepSize: n,
                    evt: t
                })
            },
            _onmousein: function(t) {
                var e = this.g.zoomer.get("stepSize" * t.rowPos),
                    n = this.g.zoomer.get("stepSize");
                return this.g.trigger("column:mousein", {
                    rowPos: e,
                    stepSize: n,
                    evt: t
                })
            },
            _onmouseout: function(t) {
                var e = this.g.zoomer.get("stepSize" * t.rowPos),
                    n = this.g.zoomer.get("stepSize");
                return this.g.trigger("column:mouseout", {
                    rowPos: e,
                    stepSize: n,
                    evt: t
                })
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(140),
        o = r(i),
        s = n(136),
        u = r(s),
        a = n(142),
        l = r(a),
        c = n(137),
        f = r(c),
        h = n(4),
        d = h.extend({
            initialize: function(t) {
                return this.g = t.g, this.blockEvents = !1, this.listenTo(this.g.vis, "change:header", function() {
                    return this.draw(), this.render()
                }), this.listenTo(this.g.vis, "change", this._setSpacer), this.listenTo(this.g.zoomer, "change:alignmentWidth", this._setWidth), this.listenTo(this.g.zoomer, "change:_alignmentScrollLeft", this._adjustScrollingLeft), this.listenTo(this.g.columns, "change:hidden", function() {
                    return this.draw(), this.render()
                }), this.draw(), this.g.vis.once("change:loaded", this._adjustScrollingLeft, this)
            },
            events: {
                scroll: "_sendScrollEvent"
            },
            draw: function() {
                if (this.removeViews(), this.g.vis.get("conserv")) {
                    var t = new u["default"]({
                        model: this.model,
                        g: this.g
                    });
                    t.ordering = -20, this.addView("conserv", t)
                }
                if (this.g.vis.get("markers")) {
                    var e = new o["default"]({
                        model: this.model,
                        g: this.g
                    });
                    e.ordering = -10, this.addView("marker", e)
                }
                if (this.g.vis.get("seqlogo")) {
                    var n = new l["default"]({
                        model: this.model,
                        g: this.g
                    });
                    n.ordering = -30, this.addView("seqlogo", n)
                }
                if (this.g.vis.get("gapHeader")) {
                    var r = new f["default"]({
                        model: this.model,
                        g: this.g
                    });
                    return r.ordering = -25, this.addView("gapview", r)
                }
            },
            render: function() {
                return this.renderSubviews(), this._setSpacer(), this.el.className = "biojs_msa_rheader", this.el.style.overflowX = "auto", this.el.style.display = "inline-block", this._setWidth(), this._adjustScrollingLeft(), this
            },
            _sendScrollEvent: function() {
                return this.blockEvents || this.g.zoomer.set("_alignmentScrollLeft", this.el.scrollLeft, {
                    origin: "header"
                }), this.blockEvents = !1
            },
            _adjustScrollingLeft: function(t, e, n) {
                if (null == ("undefined" != typeof n && null !== n ? n.origin : void 0) || "header" !== n.origin) {
                    var r = this.g.zoomer.get("_alignmentScrollLeft");
                    return this.blockEvents = !0, this.el.scrollLeft = r
                }
            },
            _setSpacer: function() {
                return this.el.style.marginLeft = this._getLabelWidth() + "px"
            },
            _getLabelWidth: function() {
                var t = 0;
                return this.g.vis.get("leftHeader") || (t += this.g.zoomer.getLeftBlockWidth()), t
            },
            _setWidth: function() {
                return this.el.style.width = this.g.zoomer.getAlignmentWidth() + "px"
            }
        });
    e["default"] = d
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(74),
        i = n(2),
        o = n(7),
        s = i.extend({
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.zoomer, "change:alignmentWidth", this.render), this.listenTo(this.g.colorscheme, "change", function() {
                    var t = this.g.colorscheme.getSelectedScheme();
                    return this.seqlogo.changeColors(t), this.render()
                }), this.listenTo(this.g.zoomer, "change:columnWidth", function() {
                    return this.seqlogo.column_width = this.g.zoomer.get("columnWidth"), this.render()
                }), this.listenTo(this.g.stats, "reset", function() {
                    return this.draw(), this.render()
                }), this.draw()
            },
            draw: function() {
                o.removeAllChilds(this.el);
                var t = this.g.stats.conservResidue({
                    scaled: !0
                });
                t = _.map(t, function(t) {
                    return _.pick(t, function(t, e) {
                        return "-" !== e
                    })
                });
                var e = {
                        alphabet: "aa",
                        heightArr: t
                    },
                    n = this.g.colorscheme.getSelectedScheme();
                return this.seqlogo = new r({
                    model: this.model,
                    g: this.g,
                    data: e,
                    yaxis: !1,
                    scroller: !1,
                    xaxis: !1,
                    height: 100,
                    column_width: this.g.zoomer.get("columnWidth"),
                    positionMarker: !1,
                    zoom: 1,
                    el: this.el,
                    colors: n
                })
            },
            render: function() {
                return this.seqlogo.render()
            }
        });
    e["default"] = s
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(144),
        o = r(i),
        s = n(4),
        u = s.extend({
            initialize: function(t) {
                var e = this;
                return this.g = t.g, this.draw(), this.listenTo(this.g.zoomer, "change:_alignmentScrollTop", this._adjustScrollingTop), this.g.vis.once("change:loaded", this._adjustScrollingTop, this), this.listenTo(this.g.zoomer, "change:alignmentHeight", this._setHeight), this.listenTo(this.model, "change:reference", this.draw), this.listenTo(this.model, "reset add remove", function() {
                    return e.draw(), e.render()
                })
            },
            draw: function() {
                this.removeViews();
                for (var t = 0; t < this.model.length; t++)
                    if (!this.model.at(t).get("hidden")) {
                        var e = new o["default"]({
                            model: this.model.at(t),
                            g: this.g
                        });
                        e.ordering = t, this.addView("row_" + t, e)
                    }
            },
            events: {
                scroll: "_sendScrollEvent"
            },
            _sendScrollEvent: function() {
                return this.g.zoomer.set("_alignmentScrollTop", this.el.scrollTop, {
                    origin: "label"
                })
            },
            _adjustScrollingTop: function() {
                return this.el.scrollTop = this.g.zoomer.get("_alignmentScrollTop")
            },
            render: function() {
                return this.renderSubviews(), this.el.className = "biojs_msa_labelblock", this.el.style.display = "inline-block", this.el.style.verticalAlign = "top", this.el.style.overflowY = "auto", this.el.style.overflowX = "hidden", this.el.style.fontSize = this.g.zoomer.get("labelFontsize") + "px", this.el.style.lineHeight = "" + this.g.zoomer.get("labelLineHeight"), this._setHeight(), this
            },
            _setHeight: function() {
                return this.el.style.height = this.g.zoomer.get("alignmentHeight") + "px"
            }
        });
    e["default"] = u
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(145),
        o = r(i),
        s = n(146),
        u = r(s),
        a = n(4),
        l = a.extend({
            initialize: function(t) {
                return this.g = t.g, this.draw(), this.listenTo(this.g.vis, "change:labels", this.drawR), this.listenTo(this.g.vis, "change:metacell", this.drawR), this.listenTo(this.g.zoomer, "change:rowHeight", function() {
                    return this.el.style.height = this.g.zoomer.get("rowHeight") + "px"
                }), this.listenTo(this.g.selcol, "change reset add", this.setSelection)
            },
            draw: function() {
                if (this.removeViews(), this.g.vis.get("labels") && this.addView("labels", new o["default"]({
                        model: this.model,
                        g: this.g
                    })), this.g.vis.get("metacell")) {
                    var t = new u["default"]({
                        model: this.model,
                        g: this.g
                    });
                    return this.addView("metacell", t)
                }
            },
            drawR: function() {
                return this.draw(), this.render()
            },
            render: function() {
                return this.renderSubviews(), this.el.setAttribute("class", "biojs_msa_labelrow"), this.el.style.height = this.g.zoomer.get("rowHeight") * (this.model.attributes.height || 1) + "px", this.setSelection(), this
            },
            setSelection: function() {
                return this.g.selcol.getSelForRow(this.model.id).length > 0 ? this.el.style.fontWeight = "bold" : this.el.style.fontWeight = "normal"
            }
        });
    e["default"] = l
}, function(t, e, n) {
    "use strict";
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var r = n(2),
        i = n(7),
        o = r.extend({
            initialize: function(t) {
                return this.seq = t.seq, this.g = t.g, this.manageEvents()
            },
            manageEvents: function() {
                var t = {};
                return this.g.config.get("registerMouseClicks") && (t.click = "_onclick"), this.g.config.get("registerMouseHover") && (t.mousein = "_onmousein", t.mouseout = "_onmouseout"), this.delegateEvents(t), this.listenTo(this.g.config, "change:registerMouseHover", this.manageEvents), this.listenTo(this.g.config, "change:registerMouseClick", this.manageEvents), this.listenTo(this.g.vis, "change:labelName change:labelId change:labelPartition change:labelCheckbox", this.render), this.listenTo(this.g.zoomer, "change:labelIdLength change:labelNameLength change:labelPartLength change:labelCheckLength", this.render), this.listenTo(this.g.zoomer, "change:labelFontSize change:labelLineHeight change:labelWidth change:rowHeight", this.render)
            },
            render: function() {
                if (i.removeAllChilds(this.el), this.el.style.width = this.g.zoomer.getLabelWidth() + "px", this.el.setAttribute("class", "biojs_msa_labels"), this.g.vis.get("labelCheckbox")) {
                    var t = document.createElement("input");
                    t.setAttribute("type", "checkbox"), t.value = this.model.get("id"), t.name = "seq", t.style.width = this.g.zoomer.get("labelCheckLength") + "px", this.el.appendChild(t)
                }
                if (this.g.vis.get("labelId")) {
                    var e = document.createElement("span"),
                        n = this.model.get("id");
                    isNaN(n) || n++, e.textContent = n, e.style.width = this.g.zoomer.get("labelIdLength") + "px", e.style.display = "inline-block", this.el.appendChild(e)
                }
                if (this.g.vis.get("labelPartition")) {
                    var r = document.createElement("span");
                    r.style.width = this.g.zoomer.get("labelPartLength") + "px", r.textContent = this.model.get("partition"), r.style.display = "inline-block", this.el.appendChild(e), this.el.appendChild(r)
                }
                if (this.g.vis.get("labelName")) {
                    var o = document.createElement("span");
                    o.textContent = this.model.get("name"), this.model.get("ref") && this.g.config.get("hasRef") && (o.style.fontWeight = "bold"), o.style.width = this.g.zoomer.get("labelNameLength") + "px", this.el.appendChild(o)
                }
                return this.el.style.overflow = scroll, this.el.style.fontSize = this.g.zoomer.get("labelFontsize") + "px", this
            },
            _onclick: function(t) {
                var e = this.model.get("id");
                return this.g.trigger("row:click", {
                    seqId: e,
                    evt: t
                })
            },
            _onmousein: function(t) {
                var e = this.model.get("id");
                return this.g.trigger("row:mouseout", {
                    seqId: e,
                    evt: t
                })
            },
            _onmouseout: function(t) {
                var e = this.model.get("id");
                return this.g.trigger("row:mouseout", {
                    seqId: e,
                    evt: t
                })
            }
        });
    e["default"] = o
}, function(t, e, n) {
    "use strict";

    function r(t) {
        return t && t.__esModule ? t : {
            "default": t
        }
    }
    Object.defineProperty(e, "__esModule", {
        value: !0
    });
    var i = n(10),
        o = n(6),
        s = r(o),
        u = n(3),
        a = n(2),
        l = n(7),
        c = a.extend({
            className: "biojs_msa_metaview",
            initialize: function(t) {
                return this.g = t.g, this.listenTo(this.g.vis, "change:metacell", this.render), this.listenTo(this.g.zoomer, "change:metaWidth", this.render)
            },
            events: {
                click: "_onclick",
                mousein: "_onmousein",
                mouseout: "_onmouseout"
            },
            render: function() {
                l.removeAllChilds(this.el), this.el.style.display = "inline-block";
                var t = this.g.zoomer.getMetaWidth();
                if (this.el.style.width = t - 10, this.el.style.paddingRight = 5, this.el.style.paddingLeft = 5, this.el.style.fontSize = this.g.zoomer.get("labelFontsize") - 2 + "px", this.g.vis.get("metaGaps")) {
                    var e = this.model.get("seq"),
                        n = (0, u.reduce)(e, function(t, e) {
                            return "-" === e ? ++t : void 0;
                        }, 0);
                    n = (100 * n / e.length).toFixed(0) + "%";
                    var r = document.createElement("span");
                    r.textContent = n, r.style.display = "inline-block", r.style.width = 35, this.el.appendChild(r)
                }
                if (this.g.vis.get("metaIdentity")) {
                    var o = this.g.stats.identity()[this.model.id],
                        a = document.createElement("span");
                    this.model.get("ref") && this.g.config.get("hasRef") ? a.textContent = "ref." : "undefined" != typeof o && null !== o && (a.textContent = o.toFixed(2)), a.style.display = "inline-block", a.style.width = 40, this.el.appendChild(a)
                }
                if (this.g.vis.get("metaLinks") && this.model.attributes.ids) {
                    var c = i.seqs.buildLinks(this.model.attributes.ids);
                    if (Object.keys(c).length > 0) {
                        var f = new s["default"]({
                            name: "↗"
                        });
                        c.forEach(function(t, e) {
                            return f.addNode(e, function(e) {
                                return window.open(t)
                            })
                        });
                        var h = f.buildDOM();
                        return h.style.cursor = "pointer", this.el.appendChild(h)
                    }
                }
            },
            _onclick: function(t) {
                return this.g.trigger("meta:click", {
                    seqId: this.model.get("id", {
                        evt: t
                    })
                })
            },
            _onmousein: function(t) {
                return this.g.trigger("meta:mousein", {
                    seqId: this.model.get("id", {
                        evt: t
                    })
                })
            },
            _onmouseout: function(t) {
                return this.g.trigger("meta:mouseout", {
                    seqId: this.model.get("id", {
                        evt: t
                    })
                })
            }
        });
    e["default"] = c
}, function(t, e) {
    (function(e) {
        "use strict";
        "undefined" != typeof window ? t.exports = window : "undefined" != typeof e ? t.exports = e : "undefined" != typeof self ? t.exports = self : t.exports = {}
    }).call(e, function() {
        return this
    }())
}, function(t, e, n) {
    "use strict";

    function r(t, e, n) {
        if (!u(e)) throw new TypeError("iterator must be a function");
        arguments.length < 3 && (n = this), "[object Array]" === a.call(t) ? i(t, e, n) : "string" == typeof t ? o(t, e, n) : s(t, e, n)
    }

    function i(t, e, n) {
        for (var r = 0, i = t.length; r < i; r++) l.call(t, r) && e.call(n, t[r], r, t)
    }

    function o(t, e, n) {
        for (var r = 0, i = t.length; r < i; r++) e.call(n, t.charAt(r), r, t)
    }

    function s(t, e, n) {
        for (var r in t) l.call(t, r) && e.call(n, t[r], r, t)
    }
    var u = n(52);
    t.exports = r;
    var a = Object.prototype.toString,
        l = Object.prototype.hasOwnProperty
}, function(t, e) {
    "use strict";

    function n(t) {
        return t.replace(/^\s*|\s*$/g, "")
    }
    e = t.exports = n, e.left = function(t) {
        return t.replace(/^\s*/, "")
    }, e.right = function(t) {
        return t.replace(/\s*$/, "")
    }
}, function(t, e, n) {
    "use strict";
    var r = n(149),
        i = n(148),
        o = function(t) {
            return "[object Array]" === Object.prototype.toString.call(t)
        };
    t.exports = function(t) {
        if (!t) return {};
        var e = {};
        return i(r(t).split("\n"), function(t) {
            var n = t.indexOf(":"),
                i = r(t.slice(0, n)).toLowerCase(),
                s = r(t.slice(n + 1));
            "undefined" == typeof e[i] ? e[i] = s : o(e[i]) ? e[i].push(s) : e[i] = [e[i], s]
        }), e
    }
}, function(t, e) {
    "use strict";

    function n() {
        for (var t = {}, e = 0; e < arguments.length; e++) {
            var n = arguments[e];
            for (var i in n) r.call(n, i) && (t[i] = n[i])
        }
        return t
    }
    t.exports = n;
    var r = Object.prototype.hasOwnProperty
}, function(t, e, n) {
    e = t.exports = n(89)(), e.push([t.id, ".biojs_msa_stage{cursor:default;line-height:normal;font-family:Helvetica}.biojs_msa_seqblock{cursor:move}.biojs_msa_layer{display:block;white-space:nowrap}.biojs_msa_labels{color:#000;display:inline-block;cursor:pointer;vertical-align:middle;overflow:hidden;text-overflow:clip}.biojs_msa_header,.biojs_msa_labels{white-space:nowrap;text-align:left}.biojs_msa_labelrow:before{content:'';display:inline-block;width:0;height:100%;vertical-align:middle}.biojs_msa_labelrow{height:100%}.biojs_msa_labelblock::-webkit-scrollbar,.biojs_msa_rheader::-webkit-scrollbar{//:none;width:7px;height:7px}.biojs_msa_labelblock::-webkit-scrollbar-thumb,.biojs_msa_rheader::-webkit-scrollbar-thumb{border-radius:4px;background-color:rgba(0,0,0,.5);box-shadow:0 0 1px hsla(0,0%,100%,.5)}.biojs_msa_marker{color:#999;white-space:nowrap}.biojs_msa_marker .msa-col-header{cursor:pointer;text-align:center}.biojs_msa_marker .msa-col-header:hover{color:red}.smenubar .smenubar_alink{background:#3498db;background-image:-webkit-linear-gradient(top,#3498db,#2980b9);background-image:linear-gradient(180deg,#3498db,#2980b9);border-radius:28px;font-family:Arial;color:#fff;padding:3px 10px;margin-left:10px;text-decoration:none}.smenubar{display:inline-block}.smenubar .smenubar_alink:hover{cursor:pointer}.smenu-dropdown{position:absolute;z-index:9999999;display:none}.smenu-dropdown .smenu-dropdown-menu,.smenu-dropdown .smenu-dropdown-panel{min-width:160px;max-width:360px;list-style:none;background:#fff;border:1px solid #ddd;border:1px solid rgba(0,0,0,.2);border-radius:6px;box-shadow:0 5px 10px rgba(0,0,0,.2);overflow:visible;padding:4px 0;margin:0}.smenu-dropdown .smenu-dropdown-panel{padding:10px}.smenu-dropdown.smenu-dropdown-scroll .smenu-dropdown-menu,.smenu-dropdown.smenu-dropdown-scroll .smenu-dropdown-panel{max-height:358px;overflow:auto}.smenu-dropdown .smenu-dropdown-menu LI{list-style:none;padding:0;margin:0;line-height:18px}.smenu-dropdown .smenu-dropdown-menu LABEL,.smenu-dropdown .smenu-dropdown-menu LI{display:block;color:#555;text-decoration:none;line-height:18px;padding:3px 15px;white-space:nowrap}.smenu-dropdown .smenu-dropdown-menu LABEL:hover,.smenu-dropdown .smenu-dropdown-menu LI:hover{background-color:#08c;color:#fff;cursor:pointer}.smenu-dropdown .smenu-dropdown-menu .smenu-dropdown-divider{font-size:1px;border-top:1px solid #e5e5e5;padding:0;margin:5px 0}.biojs_msa_div{position:relative}.biojs_msa_scale{position:absolute;bottom:0;right:0;background-color:#fff;box-shadow:0 2px 3px #999;border-radius:3px;margin:5px 0 0 auto;padding:5px;text-align:center}.biojs_msa_scale .msa-btngroup{margin:5px auto 0}.biojs_msa_scale [type=range]{cursor:pointer}.biojs_msa_scale .msa-btn-close{text-align:right;font-size:.8em;padding:2px}.biojs_msa_scale .msa-btn-open{background-color:#fff}.biojs_msa_scale .msa-hide{display:none}.msa-btn{cursor:pointer;font-size:1.1em;display:inline-block;padding:2px 8px;margin-bottom:0;border:1px solid transparent;border-radius:4px;box-sizing:border-box}.msa-btn:hover{background-color:#ddd}", ""])
}, function(t, e, n) {
    function r(t, e) {
        for (var n = 0; n < t.length; n++) {
            var r = t[n],
                i = d[r.id];
            if (i) {
                i.refs++;
                for (var o = 0; o < i.parts.length; o++) i.parts[o](r.parts[o]);
                for (; o < r.parts.length; o++) i.parts.push(l(r.parts[o], e))
            } else {
                for (var s = [], o = 0; o < r.parts.length; o++) s.push(l(r.parts[o], e));
                d[r.id] = {
                    id: r.id,
                    refs: 1,
                    parts: s
                }
            }
        }
    }

    function i(t) {
        for (var e = [], n = {}, r = 0; r < t.length; r++) {
            var i = t[r],
                o = i[0],
                s = i[1],
                u = i[2],
                a = i[3],
                l = {
                    css: s,
                    media: u,
                    sourceMap: a
                };
            n[o] ? n[o].parts.push(l) : e.push(n[o] = {
                id: o,
                parts: [l]
            })
        }
        return e
    }

    function o(t, e) {
        var n = v(),
            r = _[_.length - 1];
        if ("top" === t.insertAt) r ? r.nextSibling ? n.insertBefore(e, r.nextSibling) : n.appendChild(e) : n.insertBefore(e, n.firstChild), _.push(e);
        else {
            if ("bottom" !== t.insertAt) throw new Error("Invalid value for parameter 'insertAt'. Must be 'top' or 'bottom'.");
            n.appendChild(e)
        }
    }

    function s(t) {
        t.parentNode.removeChild(t);
        var e = _.indexOf(t);
        e >= 0 && _.splice(e, 1)
    }

    function u(t) {
        var e = document.createElement("style");
        return e.type = "text/css", o(t, e), e
    }

    function a(t) {
        var e = document.createElement("link");
        return e.rel = "stylesheet", o(t, e), e
    }

    function l(t, e) {
        var n, r, i;
        if (e.singleton) {
            var o = y++;
            n = m || (m = u(e)), r = c.bind(null, n, o, !1), i = c.bind(null, n, o, !0)
        } else t.sourceMap && "function" == typeof URL && "function" == typeof URL.createObjectURL && "function" == typeof URL.revokeObjectURL && "function" == typeof Blob && "function" == typeof btoa ? (n = a(e), r = h.bind(null, n), i = function() {
            s(n), n.href && URL.revokeObjectURL(n.href)
        }) : (n = u(e), r = f.bind(null, n), i = function() {
            s(n)
        });
        return r(t),
            function(e) {
                if (e) {
                    if (e.css === t.css && e.media === t.media && e.sourceMap === t.sourceMap) return;
                    r(t = e)
                } else i()
            }
    }

    function c(t, e, n, r) {
        var i = n ? "" : r.css;
        if (t.styleSheet) t.styleSheet.cssText = b(e, i);
        else {
            var o = document.createTextNode(i),
                s = t.childNodes;
            s[e] && t.removeChild(s[e]), s.length ? t.insertBefore(o, s[e]) : t.appendChild(o)
        }
    }

    function f(t, e) {
        var n = e.css,
            r = e.media;
        if (r && t.setAttribute("media", r), t.styleSheet) t.styleSheet.cssText = n;
        else {
            for (; t.firstChild;) t.removeChild(t.firstChild);
            t.appendChild(document.createTextNode(n))
        }
    }

    function h(t, e) {
        var n = e.css,
            r = e.sourceMap;
        r && (n += "\n/*# sourceMappingURL=data:application/json;base64," + btoa(unescape(encodeURIComponent(JSON.stringify(r)))) + " */");
        var i = new Blob([n], {
                type: "text/css"
            }),
            o = t.href;
        t.href = URL.createObjectURL(i), o && URL.revokeObjectURL(o)
    }
    var d = {},
        p = function(t) {
            var e;
            return function() {
                return "undefined" == typeof e && (e = t.apply(this, arguments)), e
            }
        },
        g = p(function() {
            return /msie [6-9]\b/.test(window.navigator.userAgent.toLowerCase())
        }),
        v = p(function() {
            return document.head || document.getElementsByTagName("head")[0]
        }),
        m = null,
        y = 0,
        _ = [];
    t.exports = function(t, e) {
        e = e || {}, "undefined" == typeof e.singleton && (e.singleton = g()), "undefined" == typeof e.insertAt && (e.insertAt = "bottom");
        var n = i(t);
        return r(n, e),
            function(t) {
                for (var o = [], s = 0; s < n.length; s++) {
                    var u = n[s],
                        a = d[u.id];
                    a.refs--, o.push(a)
                }
                t && r(i(t), e);
                for (var s = 0; s < o.length; s++) {
                    var a = o[s];
                    if (0 === a.refs) {
                        for (var l = 0; l < a.parts.length; l++) a.parts[l]();
                        delete d[a.id]
                    }
                }
            }
    };
    var b = function() {
        var t = [];
        return function(e, n) {
            return t[e] = n, t.filter(Boolean).join("\n")
        }
    }()
}, function(t, e, n) {
    var r = n(152);
    "string" == typeof r && (r = [
        [t.id, r, ""]
    ]), n(153)(r, {}), r.locals && (t.exports = r.locals)
}, function(t, e, n) {
    function r(t) {
        return n(i(t))
    }

    function i(t) {
        return o[t] || function() {
            throw new Error("Cannot find module '" + t + "'.")
        }()
    }
    var o = {
        "./StageScale": 25,
        "./StageScale.js": 25,
        "./colorscheme": 26,
        "./colorscheme.js": 26,
        "./columns": 27,
        "./columns.js": 27,
        "./config": 28,
        "./config.js": 28,
        "./package": 29,
        "./package.js": 29,
        "./selection/Selection": 8,
        "./selection/Selection.js": 8,
        "./selection/SelectionCol": 9,
        "./selection/SelectionCol.js": 9,
        "./selection/index": 49,
        "./selection/index.js": 49,
        "./user": 30,
        "./user.js": 30,
        "./visOrdering": 31,
        "./visOrdering.js": 31,
        "./visibility": 32,
        "./visibility.js": 32,
        "./zoomer": 33,
        "./zoomer.js": 33
    };
    r.keys = function() {
        return Object.keys(o)
    }, r.resolve = i, t.exports = r, r.id = 155
}, function(t, e, n) {
    function r(t) {
        return n(i(t))
    }

    function i(t) {
        return o[t] || function() {
            throw new Error("Cannot find module '" + t + "'.")
        }()
    }
    var o = {
        "./bmath": 15,
        "./bmath.js": 15,
        "./exporter": 16,
        "./exporter.js": 16,
        "./file": 17,
        "./file.js": 17,
        "./index": 37,
        "./index.js": 37,
        "./loader": 18,
        "./loader.js": 18,
        "./proxy": 19,
        "./proxy.js": 19,
        "./recognize": 38,
        "./recognize.js": 38,
        "./seqgen": 39,
        "./seqgen.js": 39,
        "./svg": 12,
        "./svg.js": 12,
        "./tree": 40,
        "./tree.js": 40
    };
    r.keys = function() {
        return Object.keys(o)
    }, r.resolve = i, t.exports = r, r.id = 156
}]);
//# sourceMappingURL=msa.js.map