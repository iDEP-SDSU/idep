(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
})(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

//ShinyGO
//ga('create', 'UA-87863704-1', 'auto');
//iDEP at DO
ga('create', 'UA-87863704-2', 'auto');
//iDEP at prod
//ga('create', 'UA-87863704-3', 'auto');
ga('send', 'pageview');
// Event Tracking Code
$(document).on('shiny:value', function(event) {
//    if(event.name == 'navBar'){
        ga('send', 'pageview');
//    }
});