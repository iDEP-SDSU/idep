  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-87863704-1', 'auto');
  ga('send', 'pageview');
$(document).on('shiny:inputchanged', function(event) {
        ga('send', 'pageview');
});

$(document).on('shiny:inputchanged', function(event) {
if(event.name == 'file1'){
        ga('send', 'event','upload','expression','data');
   }
});

$(document).on('shiny:error', function(event) {
        ga('send', 'event','error','error','error');
});    
# dtest   