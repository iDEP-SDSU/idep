import Firebase from 'firebase'
console.log(Firebase)
let config = {
    apiKey: "AIzaSyBO4CCJzL7U9pFSEv-9ETqVt5dzMNKiwk4",
    authDomain: "bcloud.firebaseapp.com",
    databaseURL: "https://bcloud.firebaseio.com",
    projectId: "firebase-bcloud",
    storageBucket: "firebase-bcloud.appspot.com",
    messagingSenderId: "172712893865"
};

if (!Firebase.apps.length) {
    Firebase.initializeApp(config)
}

// Export the database for components to use.
// If you want to get fancy, use mixins or provide / inject to avoid redundant imports.
export default Firebase.database();
