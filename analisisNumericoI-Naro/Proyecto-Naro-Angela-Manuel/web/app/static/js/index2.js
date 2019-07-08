var container = document.getElementById('container');
setTimeout(function() {
	container.classList.add('close');
  //document.body.style.overflowY= "visible";// despueés de cargar le devolvemos el scroll
}, 4600);
var initial = document.getElementById('p1');
var button = document.getElementById('button');
setTimeout(function() {
  var name = prompt("Ingresa tu nombre.");
  var string = " "+name+"!";
  var node=document.createTextNode(string);
  initial.appendChild(node);
  button.classList.add('show');
}, 5000);