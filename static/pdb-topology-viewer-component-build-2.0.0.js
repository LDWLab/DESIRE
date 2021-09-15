!function(t){var e={};function n(r){if(e[r])return e[r].exports;var i=e[r]={i:r,l:!1,exports:{}};return t[r].call(i.exports,i,i.exports,n),i.l=!0,i.exports}n.m=t,n.c=e,n.d=function(t,e,r){n.o(t,e)||Object.defineProperty(t,e,{enumerable:!0,get:r})},n.r=function(t){"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(t,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(t,"__esModule",{value:!0})},n.t=function(t,e){if(1&e&&(t=n(t)),8&e)return t;if(4&e&&"object"==typeof t&&t&&t.__esModule)return t;var r=Object.create(null);if(n.r(r),Object.defineProperty(r,"default",{enumerable:!0,value:t}),2&e&&"string"!=typeof t)for(var i in t)n.d(r,i,function(e){return t[e]}.bind(null,i));return r},n.n=function(t){var e=t&&t.__esModule?function(){return t.default}:function(){return t};return n.d(e,"a",e),e},n.o=function(t,e){return Object.prototype.hasOwnProperty.call(t,e)},n.p="",n(n.s=7)}([function(t,e){function n(e,r){return t.exports=n=Object.setPrototypeOf||function(t,e){return t.__proto__=e,t},n(e,r)}t.exports=n},function(t,e){function n(e){return t.exports=n=Object.setPrototypeOf?Object.getPrototypeOf:function(t){return t.__proto__||Object.getPrototypeOf(t)},n(e)}t.exports=n},function(t,e){function n(t,e){for(var n=0;n<e.length;n++){var r=e[n];r.enumerable=r.enumerable||!1,r.configurable=!0,"value"in r&&(r.writable=!0),Object.defineProperty(t,r.key,r)}}t.exports=function(t,e,r){return e&&n(t.prototype,e),r&&n(t,r),t}},function(t,e){t.exports=function(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")}},function(t,e,n){var r=n(8),i=n(9);t.exports=function(t,e){return!e||"object"!==r(e)&&"function"!=typeof e?i(t):e}},function(t,e,n){var r=n(0);t.exports=function(t,e){if("function"!=typeof e&&null!==e)throw new TypeError("Super expression must either be null or a function");t.prototype=Object.create(e&&e.prototype,{constructor:{value:t,writable:!0,configurable:!0}}),e&&r(t,e)}},function(t,e,n){var r=n(1),i=n(0),o=n(10),u=n(11);function s(e){var n="function"==typeof Map?new Map:void 0;return t.exports=s=function(t){if(null===t||!o(t))return t;if("function"!=typeof t)throw new TypeError("Super expression must either be null or a function");if(void 0!==n){if(n.has(t))return n.get(t);n.set(t,e)}function e(){return u(t,arguments,r(this).constructor)}return e.prototype=Object.create(t.prototype,{constructor:{value:e,enumerable:!1,writable:!0,configurable:!0}}),i(e,t)},s(e)}t.exports=s},function(t,e,n){"use strict";n.r(e);var r=n(3),i=n.n(r),o=n(4),u=n.n(o),s=n(1),l=n.n(s),c=n(2),f=n.n(c),a=n(5),p=n.n(a),y=n(6),h=function(t){function e(){return i()(this,e),u()(this,l()(e).call(this))}return p()(e,t),f()(e,null,[{key:"observedAttributes",get:function(){return["entry-id","entity-id","filter-range","chain-id","display-style","error-style","menu-style","subscribe-events","pvapi"]}}]),f()(e,[{key:"validateParams",value:function(){return void 0!==this.entryId&&void 0!==this.entityId&&null!=this.entryId&&null!=this.entityId}},{key:"invokePlugin",value:function(){if(this.validateParams()){void 0===this.pluginInstance&&(this.pluginInstance=new PdbTopologyViewerPlugin);var t={entryId:this.entryId,entityId:this.entityId,entropyId:this.entropyId};void 0!==this.chainId&&null!==this.chainId&&(t.chainId=this.chainId),void 0!==this.displayStyle&&null!==this.displayStyle&&(t.displayStyle=this.displayStyle),void 0!==this.errorStyle&&null!==this.errorStyle&&(t.errorStyle=this.errorStyle),void 0!==this.menuStyle&&null!==this.menuStyle&&(t.menuStyle=this.menuStyle),void 0!==this.subscribeEvents&&null!==this.subscribeEvents&&(t.subscribeEvents=this.subscribeEvents),void 0!==this.pvAPI&&null!==this.pvAPI&&(t.pvAPI=this.pvAPI),void 0!==this.filterRange&&null!==this.filterRange&&(t.filterRange=this.filterRange),this.pluginInstance.render(this,t)}}},{key:"attributeChangedCallback",value:function(){this.entryId=this.getAttribute("entry-id"),this.entityId=this.getAttribute("entity-id"),this.chainId=this.getAttribute("chain-id"),this.filterRange=this.getAttribute("filter-range"),this.displayStyle=this.getAttribute("display-style"),this.errorStyle=this.getAttribute("error-style"),this.menuStyle=this.getAttribute("menu-style"),this.subscribeEvents=this.getAttribute("subscribe-events"),this.pvAPI=/true/i.test(this.getAttribute("pvapi")),this.invokePlugin()}}]),e}(n.n(y)()(HTMLElement));e.default=h,customElements.define("pdb-topology-viewer",h)},function(t,e){function n(e){return"function"==typeof Symbol&&"symbol"==typeof Symbol.iterator?t.exports=n=function(t){return typeof t}:t.exports=n=function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":typeof t},n(e)}t.exports=n},function(t,e){t.exports=function(t){if(void 0===t)throw new ReferenceError("this hasn't been initialised - super() hasn't been called");return t}},function(t,e){t.exports=function(t){return-1!==Function.toString.call(t).indexOf("[native code]")}},function(t,e,n){var r=n(0);function i(){if("undefined"==typeof Reflect||!Reflect.construct)return!1;if(Reflect.construct.sham)return!1;if("function"==typeof Proxy)return!0;try{return Date.prototype.toString.call(Reflect.construct(Date,[],(function(){}))),!0}catch(t){return!1}}function o(e,n,u){return i()?t.exports=o=Reflect.construct:t.exports=o=function(t,e,n){var i=[null];i.push.apply(i,e);var o=new(Function.bind.apply(t,i));return n&&r(o,n.prototype),o},o.apply(null,arguments)}t.exports=o}]);
//# sourceMappingURL=pdb-topology-viewer-component-build-2.0.0.js.map