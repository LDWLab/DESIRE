import Vue from 'vue';
import RV3 from './RV3.vue';
import VueTour from 'vue-tour'

function loadScript(src) {
  return new Promise(function (resolve, reject) {
      var s;
      s = document.createElement('script');
      s.src = src;
      s.onload = resolve;
      s.onerror = reject;
      document.head.appendChild(s);
  });
}
loadScript('/static/alignments/RV3_helpers.js')
//.catch(loadScript.bind(null, localSource))
//.then(successCallback, failureCallback);
//Use these two to catch failures
loadScript('/static/alignments/RV3_after.js')

//From here https://github.com/pulsardev/vue-tour
Vue.use(VueTour)

var vm_both = new Vue({
  el: '#app',
  template: '<RV3/>',
  components: { RV3 }
})

window.vm = vm_both.$children[0].$children[1];