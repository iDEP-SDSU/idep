<template>
  <section class="container">
    <div>
      <logo/>
      <h1 class="title">
        idep-loader
      </h1>
      <h2 class="subtitle">
        Welcome to iDep sss
      </h2>
      <h2 class="subtitle">
        You will be redirected to ideal iDep application. 
        {{message}}
      </h2>
      <ul class="list">
        <li v-for="idep in ideps">
          <div class="f4">{{idep.id}}</div> <div v-for="(conn, idx) in idep.connections">{{idx}}</div>
        </li>
      </ul>
    </div>
  </section>
</template>

<script>
import Logo from '~/components/Logo.vue'
import fb from "~/plugins/fbConn"

export default {
  components: {
    Logo
  },
  data() {
    return{
      ideps: [],
      message: ""
    }
  },
  mounted : function() {
    var vm = this;
      var usersShinyRef = fb.ref("usersShiny").once("value").then(function(snap){
        vm.ideps = snap.val()
        snap.forEach(idep =>{
          let crntIdep = idep.val()
          if(crntIdep.connections == undefined){
            console.log(crntIdep.id, "free")
            vm.message = crntIdep.id + " is ideal. After 3 seconds, you will be redirected "
            setInterval(function(){ 
              window.location.href = "http://bioinformatics.sdstate.edu/"+crntIdep.id;
            }, 3000);
            
          }
          else{
            console.log(crntIdep.id, length(crntIdep.connections))
          } 
        })
        
    })      
  }
}
</script>

<style>
  .container {
    min-height: 100vh;
    display: flex;
    justify-content: center;
    align-items: center;
    text-align: center;
  }

  .title {
    font-family: "Quicksand", "Source Sans Pro", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; /* 1 */
    display: block;
    font-weight: 300;
    font-size: 100px;
    color: #35495e;
    letter-spacing: 1px;
  }

  .subtitle {
    font-weight: 300;
    font-size: 42px;
    color: #526488;
    word-spacing: 5px;
    padding-bottom: 15px;
  }

  .links {
    padding-top: 15px;
  }
</style>