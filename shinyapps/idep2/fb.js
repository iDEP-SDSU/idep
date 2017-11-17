function guid() {
  function s4() {
    return Math.floor((1 + Math.random()) * 0x10000)
      .toString(16)
      .substring(1);
  }
  return s4() + s4() + '-' + s4() + '-' + s4() + '-' +
    s4() + '-' + s4() + s4() + s4();
}
var idepName = 'idep2'
var config = {
    apiKey: "AIzaSyBO4CCJzL7U9pFSEv-9ETqVt5dzMNKiwk4",
    authDomain: "bcloud.firebaseapp.com",
    databaseURL: "https://bcloud.firebaseio.com",
    projectId: "firebase-bcloud",
    storageBucket: "firebase-bcloud.appspot.com",
    messagingSenderId: "172712893865"
    };
firebase.initializeApp(config);
var fb = firebase.database();
var uuid = guid()

var myConnectionsRef = firebase.database().ref('usersShiny/'+idepName+'/connections');
// stores the timestamp of my last disconnect (the last time I was seen online)
var lastOnlineRef = firebase.database().ref('usersShiny/'+idepName+'/lastOnline');
var infoRef = firebase.database().ref('usersShiny/'+idepName+'/id');
var connectedRef = firebase.database().ref('.info/connected');
connectedRef.on('value', function(snap) {
  console.log("this not happen")
  if (snap.val() === true) {
    infoRef.set(idepName)
    console.log("this happen")
    // We're connected (or reconnected)! Do anything here that should happen only if online (or on reconnect)
    var con = myConnectionsRef.push();
    // When I disconnect, remove this device
    con.onDisconnect().remove();
    // Add this device to my connections list
    // this value could contain info about the device or a timestamp too
    con.set(true);
    // When I disconnect, update the last time I was seen online
    lastOnlineRef.onDisconnect().set(firebase.database.ServerValue.TIMESTAMP);
  }
});
var orgRef = fb.ref("idep/gmt/org")
orgRef.on('value', function(snap){
  console.log(snap.val())
})