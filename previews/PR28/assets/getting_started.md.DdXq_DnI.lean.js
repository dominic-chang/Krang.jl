import{_ as n}from"./chunks/coordinate.D_JvmGHt.js";import{_ as o}from"./chunks/raytrace.BN_57Oz_.js";import{_ as r,c as s,a2 as l,j as t,a,o as i}from"./chunks/framework.DHkBfLL0.js";const D=JSON.parse('{"title":"What is Krang.jl?","description":"","frontmatter":{},"headers":[],"relativePath":"getting_started.md","filePath":"getting_started.md","lastUpdated":null}'),d={name:"getting_started.md"},Q={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},T={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.355ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.999ex",height:"1.95ex",role:"img",focusable:"false",viewBox:"0 -705 883.6 862.1","aria-hidden":"true"},m={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},h={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.029ex"},xmlns:"http://www.w3.org/2000/svg",width:"3.054ex",height:"1.649ex",role:"img",focusable:"false",viewBox:"0 -716 1350 729","aria-hidden":"true"},p={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.023ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.061ex",height:"1.618ex",role:"img",focusable:"false",viewBox:"0 -705 469 715","aria-hidden":"true"},c={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},u={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.029ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.17ex",height:"1.005ex",role:"img",focusable:"false",viewBox:"0 -431 517 444","aria-hidden":"true"},k={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},x={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.029ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.92ex",height:"1.747ex",role:"img",focusable:"false",viewBox:"0 -759 848.5 772","aria-hidden":"true"},y={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},w={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"15.569ex",height:"2.396ex",role:"img",focusable:"false",viewBox:"0 -809 6881.7 1059","aria-hidden":"true"},f={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},b={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"22.916ex",height:"2.565ex",role:"img",focusable:"false",viewBox:"0 -883.9 10128.7 1133.9","aria-hidden":"true"};function M(v,e,L,H,E,Z){return i(),s("div",null,[e[22]||(e[22]=l('<h1 id="What-is-Krang.jl?" tabindex="-1">What is Krang.jl? <a class="header-anchor" href="#What-is-Krang.jl?" aria-label="Permalink to &quot;What is Krang.jl? {#What-is-Krang.jl?}&quot;">​</a></h1><p>Kerr Raytracer for Analytic Null Geodesics (Krang) is a raytracing code for geometries that are embedded in the Kerr spacetime. Krang solves the Kerr geodesic problem using the analytic solutions that have been derived from the Hamilton-Jacobi Formalism. This choice makes the code efficient and accurate, and is an ideal formalism to isolated sub images that manifest from strong lensing.</p><h2 id="philosophy" tabindex="-1">Philosophy <a class="header-anchor" href="#philosophy" aria-label="Permalink to &quot;Philosophy&quot;">​</a></h2><p>The current philosophy of Krang is to raytrace meshes. Meshes are made by &#39;painting&#39; materials onto geometries. These are then viewed by cameras that are sensitive to particular observables (<em>e.g. intensity and polarization</em>).</p><h2 id="meshes" tabindex="-1">Meshes <a class="header-anchor" href="#meshes" aria-label="Permalink to &quot;Meshes&quot;">​</a></h2><h3 id="materials" tabindex="-1">Materials <a class="header-anchor" href="#materials" aria-label="Permalink to &quot;Materials&quot;">​</a></h3><p>Marterials define the local physics necessary to properly render geometries. Materials may sometimes need additional information that can be stored in geometries by passing <code>attributes</code> to the geometry constructor (see <a href="./examples/polarization-example">Custom Dual Cone Model</a> for an example of passing attributes to geometries).</p><h3 id="geometries" tabindex="-1">Geometries <a class="header-anchor" href="#geometries" aria-label="Permalink to &quot;Geometries&quot;">​</a></h3><p>There are two basic geometries currently available in Krang.</p><ul><li><p><code>ConeGeometry</code> : The first is a spin axis centered cone with its apex placed at the coordinate origin.</p></li><li><p><code>MeshGeometry</code> : Geometry made from a triangular mesh. The mesh is embdedded by placing vertices at points in a Cartesian like coordinate system generated from the Boyer-Lindquist coordinates to a Cartesian like equivalent, and . There are convenience functions defines to <code>translate</code>, <code>rotate</code> and <code>scale</code> these geometries.</p></li></ul><h2 id="raytracing" tabindex="-1">Raytracing <a class="header-anchor" href="#raytracing" aria-label="Permalink to &quot;Raytracing&quot;">​</a></h2>',11)),t("p",null,[e[4]||(e[4]=a("Light rays in this module can be parameterized in terms of either the cones (")),t("mjx-container",Q,[(i(),s("svg",T,e[0]||(e[0]=[l('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="msub"><g data-mml-node="mi"><path data-c="1D703" d="M35 200Q35 302 74 415T180 610T319 704Q320 704 327 704T339 705Q393 701 423 656Q462 596 462 495Q462 380 417 261T302 66T168 -10H161Q125 -10 99 10T60 63T41 130T35 200ZM383 566Q383 668 330 668Q294 668 260 623T204 521T170 421T157 371Q206 370 254 370L351 371Q352 372 359 404T375 484T383 566ZM113 132Q113 26 166 26Q181 26 198 36T239 74T287 161T335 307L340 324H145Q145 321 136 286T120 208T113 132Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(502,-150) scale(0.707)"><path data-c="1D460" d="M131 289Q131 321 147 354T203 415T300 442Q362 442 390 415T419 355Q419 323 402 308T364 292Q351 292 340 300T328 326Q328 342 337 354T354 372T367 378Q368 378 368 379Q368 382 361 388T336 399T297 405Q249 405 227 379T204 326Q204 301 223 291T278 274T330 259Q396 230 396 163Q396 135 385 107T352 51T289 7T195 -10Q118 -10 86 19T53 87Q53 126 74 143T118 160Q133 160 146 151T160 120Q160 94 142 76T111 58Q109 57 108 57T107 55Q108 52 115 47T146 34T201 27Q237 27 263 38T301 66T318 97T323 122Q323 150 302 164T254 181T195 196T148 231Q131 256 131 289Z" style="stroke-width:3;"></path></g></g></g></g>',1)]))),e[1]||(e[1]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("msub",null,[t("mi",null,"θ"),t("mi",null,"s")])])],-1))]),e[5]||(e[5]=a("), or the minotime (")),t("mjx-container",m,[(i(),s("svg",h,e[2]||(e[2]=[l('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="mi"><path data-c="394" d="M51 0Q46 4 46 7Q46 9 215 357T388 709Q391 716 416 716Q439 716 444 709Q447 705 616 357T786 7Q786 4 781 0H51ZM507 344L384 596L137 92L383 91H630Q630 93 507 344Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(833,0)"><path data-c="1D70F" d="M39 284Q18 284 18 294Q18 301 45 338T99 398Q134 425 164 429Q170 431 332 431Q492 431 497 429Q517 424 517 402Q517 388 508 376T485 360Q479 358 389 358T299 356Q298 355 283 274T251 109T233 20Q228 5 215 -4T186 -13Q153 -13 153 20V30L203 192Q214 228 227 272T248 336L254 357Q254 358 208 358Q206 358 197 358T183 359Q105 359 61 295Q56 287 53 286T39 284Z" style="stroke-width:3;"></path></g></g></g>',1)]))),e[3]||(e[3]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("mi",{mathvariant:"normal"},"Δ"),t("mi",null,"τ")])],-1))]),e[6]||(e[6]=a("). Parameterization in terms of cones allows for images to be divided into sub images and ray traced individually."))]),e[23]||(e[23]=t("h3",{id:"Raytracing-conical-surfaces",tabindex:"-1"},[a("Raytracing conical surfaces "),t("a",{class:"header-anchor",href:"#Raytracing-conical-surfaces","aria-label":'Permalink to "Raytracing conical surfaces {#Raytracing-conical-surfaces}"'},"​")],-1)),t("p",null,[e[9]||(e[9]=a("Surfaces of constant ")),t("mjx-container",p,[(i(),s("svg",g,e[7]||(e[7]=[t("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[t("g",{"data-mml-node":"math"},[t("g",{"data-mml-node":"mi"},[t("path",{"data-c":"1D703",d:"M35 200Q35 302 74 415T180 610T319 704Q320 704 327 704T339 705Q393 701 423 656Q462 596 462 495Q462 380 417 261T302 66T168 -10H161Q125 -10 99 10T60 63T41 130T35 200ZM383 566Q383 668 330 668Q294 668 260 623T204 521T170 421T157 371Q206 370 254 370L351 371Q352 372 359 404T375 484T383 566ZM113 132Q113 26 166 26Q181 26 198 36T239 74T287 161T335 307L340 324H145Q145 321 136 286T120 208T113 132Z",style:{"stroke-width":"3"}})])])],-1)]))),e[8]||(e[8]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("mi",null,"θ")])],-1))]),e[10]||(e[10]=a(" define spin axis centered cones whose apex lie at the origin of the Boyer-Lindquist coordinate system."))]),e[24]||(e[24]=t("p",null,[t("img",{src:n,alt:""})],-1)),e[25]||(e[25]=t("p",{style:{"text-align":"center"}},"n=0 and n=1 images of emission coordinates originating from conical surfaces.",-1)),e[26]||(e[26]=t("h3",{id:"Raytracing-with-rays-parameterized-by-Mino-time",tabindex:"-1"},[a("Raytracing with rays parameterized by Mino time "),t("a",{class:"header-anchor",href:"#Raytracing-with-rays-parameterized-by-Mino-time","aria-label":'Permalink to "Raytracing with rays parameterized by Mino time {#Raytracing-with-rays-parameterized-by-Mino-time}"'},"​")],-1)),t("p",null,[e[15]||(e[15]=a("Mino time, ")),t("mjx-container",c,[(i(),s("svg",u,e[11]||(e[11]=[t("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[t("g",{"data-mml-node":"math"},[t("g",{"data-mml-node":"mi"},[t("path",{"data-c":"1D70F",d:"M39 284Q18 284 18 294Q18 301 45 338T99 398Q134 425 164 429Q170 431 332 431Q492 431 497 429Q517 424 517 402Q517 388 508 376T485 360Q479 358 389 358T299 356Q298 355 283 274T251 109T233 20Q228 5 215 -4T186 -13Q153 -13 153 20V30L203 192Q214 228 227 272T248 336L254 357Q254 358 208 358Q206 358 197 358T183 359Q105 359 61 295Q56 287 53 286T39 284Z",style:{"stroke-width":"3"}})])])],-1)]))),e[12]||(e[12]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("mi",null,"τ")])],-1))]),e[16]||(e[16]=a(", is a parameter monotonic in affine parameter, ")),t("mjx-container",k,[(i(),s("svg",x,e[13]||(e[13]=[l('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="msup"><g data-mml-node="mi"><path data-c="1D70F" d="M39 284Q18 284 18 294Q18 301 45 338T99 398Q134 425 164 429Q170 431 332 431Q492 431 497 429Q517 424 517 402Q517 388 508 376T485 360Q479 358 389 358T299 356Q298 355 283 274T251 109T233 20Q228 5 215 -4T186 -13Q153 -13 153 20V30L203 192Q214 228 227 272T248 336L254 357Q254 358 208 358Q206 358 197 358T183 359Q105 359 61 295Q56 287 53 286T39 284Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(604,363) scale(0.707)"><path data-c="2032" d="M79 43Q73 43 52 49T30 61Q30 68 85 293T146 528Q161 560 198 560Q218 560 240 545T262 501Q262 496 260 486Q259 479 173 263T84 45T79 43Z" style="stroke-width:3;"></path></g></g></g></g>',1)]))),e[14]||(e[14]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("msup",null,[t("mi",null,"τ"),t("mo",{"data-mjx-alternate":"1"},"′")])])],-1))]),e[17]||(e[17]=a(", defined by"))]),t("mjx-container",y,[(i(),s("svg",w,e[18]||(e[18]=[l('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="mi"><path data-c="1D451" d="M366 683Q367 683 438 688T511 694Q523 694 523 686Q523 679 450 384T375 83T374 68Q374 26 402 26Q411 27 422 35Q443 55 463 131Q469 151 473 152Q475 153 483 153H487H491Q506 153 506 145Q506 140 503 129Q490 79 473 48T445 8T417 -8Q409 -10 393 -10Q359 -10 336 5T306 36L300 51Q299 52 296 50Q294 48 292 46Q233 -10 172 -10Q117 -10 75 30T33 157Q33 205 53 255T101 341Q148 398 195 420T280 442Q336 442 364 400Q369 394 369 396Q370 400 396 505T424 616Q424 629 417 632T378 637H357Q351 643 351 645T353 664Q358 683 366 683ZM352 326Q329 405 277 405Q242 405 210 374T160 293Q131 214 119 129Q119 126 119 118T118 106Q118 61 136 44T179 26Q233 26 290 98L298 109L352 326Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(520,0)"><path data-c="1D70F" d="M39 284Q18 284 18 294Q18 301 45 338T99 398Q134 425 164 429Q170 431 332 431Q492 431 497 429Q517 424 517 402Q517 388 508 376T485 360Q479 358 389 358T299 356Q298 355 283 274T251 109T233 20Q228 5 215 -4T186 -13Q153 -13 153 20V30L203 192Q214 228 227 272T248 336L254 357Q254 358 208 358Q206 358 197 358T183 359Q105 359 61 295Q56 287 53 286T39 284Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(1314.8,0)"><path data-c="3D" d="M56 347Q56 360 70 367H707Q722 359 722 347Q722 336 708 328L390 327H72Q56 332 56 347ZM56 153Q56 168 72 173H708Q722 163 722 153Q722 140 707 133H70Q56 140 56 153Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(2370.6,0)"><path data-c="3A3" d="M666 247Q664 244 652 126T638 4V0H351Q131 0 95 0T57 5V6Q54 12 57 17L73 36Q89 54 121 90T182 159L305 299L56 644L55 658Q55 677 60 681Q63 683 351 683H638V679Q640 674 652 564T666 447V443H626V447Q618 505 604 543T559 605Q529 626 478 631T333 637H294H189L293 494Q314 465 345 422Q400 346 400 340Q400 338 399 337L154 57Q407 57 428 58Q476 60 508 68T551 83T575 103Q595 125 608 162T624 225L626 251H666V247Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(3092.6,0)"><path data-c="28" d="M94 250Q94 319 104 381T127 488T164 576T202 643T244 695T277 729T302 750H315H319Q333 750 333 741Q333 738 316 720T275 667T226 581T184 443T167 250T184 58T225 -81T274 -167T316 -220T333 -241Q333 -250 318 -250H315H302L274 -226Q180 -141 137 -14T94 250Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(3481.6,0)"><path data-c="1D45F" d="M21 287Q22 290 23 295T28 317T38 348T53 381T73 411T99 433T132 442Q161 442 183 430T214 408T225 388Q227 382 228 382T236 389Q284 441 347 441H350Q398 441 422 400Q430 381 430 363Q430 333 417 315T391 292T366 288Q346 288 334 299T322 328Q322 376 378 392Q356 405 342 405Q286 405 239 331Q229 315 224 298T190 165Q156 25 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 114 189T154 366Q154 405 128 405Q107 405 92 377T68 316T57 280Q55 278 41 278H27Q21 284 21 287Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(3932.6,0)"><path data-c="2C" d="M78 35T78 60T94 103T137 121Q165 121 187 96T210 8Q210 -27 201 -60T180 -117T154 -158T130 -185T117 -194Q113 -194 104 -185T95 -172Q95 -168 106 -156T131 -126T157 -76T173 -3V9L172 8Q170 7 167 6T161 3T152 1T140 0Q113 0 96 17Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(4377.2,0)"><path data-c="1D703" d="M35 200Q35 302 74 415T180 610T319 704Q320 704 327 704T339 705Q393 701 423 656Q462 596 462 495Q462 380 417 261T302 66T168 -10H161Q125 -10 99 10T60 63T41 130T35 200ZM383 566Q383 668 330 668Q294 668 260 623T204 521T170 421T157 371Q206 370 254 370L351 371Q352 372 359 404T375 484T383 566ZM113 132Q113 26 166 26Q181 26 198 36T239 74T287 161T335 307L340 324H145Q145 321 136 286T120 208T113 132Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(4846.2,0)"><path data-c="29" d="M60 749L64 750Q69 750 74 750H86L114 726Q208 641 251 514T294 250Q294 182 284 119T261 12T224 -76T186 -143T145 -194T113 -227T90 -246Q87 -249 86 -250H74Q66 -250 63 -250T58 -247T55 -238Q56 -237 66 -225Q221 -64 221 250T66 725Q56 737 55 738Q55 746 60 749Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(5235.2,0)"><path data-c="1D451" d="M366 683Q367 683 438 688T511 694Q523 694 523 686Q523 679 450 384T375 83T374 68Q374 26 402 26Q411 27 422 35Q443 55 463 131Q469 151 473 152Q475 153 483 153H487H491Q506 153 506 145Q506 140 503 129Q490 79 473 48T445 8T417 -8Q409 -10 393 -10Q359 -10 336 5T306 36L300 51Q299 52 296 50Q294 48 292 46Q233 -10 172 -10Q117 -10 75 30T33 157Q33 205 53 255T101 341Q148 398 195 420T280 442Q336 442 364 400Q369 394 369 396Q370 400 396 505T424 616Q424 629 417 632T378 637H357Q351 643 351 645T353 664Q358 683 366 683ZM352 326Q329 405 277 405Q242 405 210 374T160 293Q131 214 119 129Q119 126 119 118T118 106Q118 61 136 44T179 26Q233 26 290 98L298 109L352 326Z" style="stroke-width:3;"></path></g><g data-mml-node="msup" transform="translate(5755.2,0)"><g data-mml-node="mi"><path data-c="1D70F" d="M39 284Q18 284 18 294Q18 301 45 338T99 398Q134 425 164 429Q170 431 332 431Q492 431 497 429Q517 424 517 402Q517 388 508 376T485 360Q479 358 389 358T299 356Q298 355 283 274T251 109T233 20Q228 5 215 -4T186 -13Q153 -13 153 20V30L203 192Q214 228 227 272T248 336L254 357Q254 358 208 358Q206 358 197 358T183 359Q105 359 61 295Q56 287 53 286T39 284Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(604,413) scale(0.707)"><path data-c="2032" d="M79 43Q73 43 52 49T30 61Q30 68 85 293T146 528Q161 560 198 560Q218 560 240 545T262 501Q262 496 260 486Q259 479 173 263T84 45T79 43Z" style="stroke-width:3;"></path></g></g><g data-mml-node="mo" transform="translate(6603.7,0)"><path data-c="2C" d="M78 35T78 60T94 103T137 121Q165 121 187 96T210 8Q210 -27 201 -60T180 -117T154 -158T130 -185T117 -194Q113 -194 104 -185T95 -172Q95 -168 106 -156T131 -126T157 -76T173 -3V9L172 8Q170 7 167 6T161 3T152 1T140 0Q113 0 96 17Z" style="stroke-width:3;"></path></g></g></g>',1)]))),e[19]||(e[19]=t("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[t("mi",null,"d"),t("mi",null,"τ"),t("mo",null,"="),t("mi",{mathvariant:"normal"},"Σ"),t("mo",{stretchy:"false"},"("),t("mi",null,"r"),t("mo",null,","),t("mi",null,"θ"),t("mo",{stretchy:"false"},")"),t("mi",null,"d"),t("msup",null,[t("mi",null,"τ"),t("mo",{"data-mjx-alternate":"1"},"′")]),t("mo",null,",")])],-1))]),e[27]||(e[27]=t("p",null,"where",-1)),t("mjx-container",f,[(i(),s("svg",b,e[20]||(e[20]=[l('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="mi"><path data-c="3A3" d="M666 247Q664 244 652 126T638 4V0H351Q131 0 95 0T57 5V6Q54 12 57 17L73 36Q89 54 121 90T182 159L305 299L56 644L55 658Q55 677 60 681Q63 683 351 683H638V679Q640 674 652 564T666 447V443H626V447Q618 505 604 543T559 605Q529 626 478 631T333 637H294H189L293 494Q314 465 345 422Q400 346 400 340Q400 338 399 337L154 57Q407 57 428 58Q476 60 508 68T551 83T575 103Q595 125 608 162T624 225L626 251H666V247Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(722,0)"><path data-c="28" d="M94 250Q94 319 104 381T127 488T164 576T202 643T244 695T277 729T302 750H315H319Q333 750 333 741Q333 738 316 720T275 667T226 581T184 443T167 250T184 58T225 -81T274 -167T316 -220T333 -241Q333 -250 318 -250H315H302L274 -226Q180 -141 137 -14T94 250Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(1111,0)"><path data-c="1D45F" d="M21 287Q22 290 23 295T28 317T38 348T53 381T73 411T99 433T132 442Q161 442 183 430T214 408T225 388Q227 382 228 382T236 389Q284 441 347 441H350Q398 441 422 400Q430 381 430 363Q430 333 417 315T391 292T366 288Q346 288 334 299T322 328Q322 376 378 392Q356 405 342 405Q286 405 239 331Q229 315 224 298T190 165Q156 25 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 114 189T154 366Q154 405 128 405Q107 405 92 377T68 316T57 280Q55 278 41 278H27Q21 284 21 287Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(1562,0)"><path data-c="2C" d="M78 35T78 60T94 103T137 121Q165 121 187 96T210 8Q210 -27 201 -60T180 -117T154 -158T130 -185T117 -194Q113 -194 104 -185T95 -172Q95 -168 106 -156T131 -126T157 -76T173 -3V9L172 8Q170 7 167 6T161 3T152 1T140 0Q113 0 96 17Z" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(2006.7,0)"><path data-c="1D703" d="M35 200Q35 302 74 415T180 610T319 704Q320 704 327 704T339 705Q393 701 423 656Q462 596 462 495Q462 380 417 261T302 66T168 -10H161Q125 -10 99 10T60 63T41 130T35 200ZM383 566Q383 668 330 668Q294 668 260 623T204 521T170 421T157 371Q206 370 254 370L351 371Q352 372 359 404T375 484T383 566ZM113 132Q113 26 166 26Q181 26 198 36T239 74T287 161T335 307L340 324H145Q145 321 136 286T120 208T113 132Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(2475.7,0)"><path data-c="29" d="M60 749L64 750Q69 750 74 750H86L114 726Q208 641 251 514T294 250Q294 182 284 119T261 12T224 -76T186 -143T145 -194T113 -227T90 -246Q87 -249 86 -250H74Q66 -250 63 -250T58 -247T55 -238Q56 -237 66 -225Q221 -64 221 250T66 725Q56 737 55 738Q55 746 60 749Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(3142.4,0)"><path data-c="3D" d="M56 347Q56 360 70 367H707Q722 359 722 347Q722 336 708 328L390 327H72Q56 332 56 347ZM56 153Q56 168 72 173H708Q722 163 722 153Q722 140 707 133H70Q56 140 56 153Z" style="stroke-width:3;"></path></g><g data-mml-node="msup" transform="translate(4198.2,0)"><g data-mml-node="mi"><path data-c="1D45F" d="M21 287Q22 290 23 295T28 317T38 348T53 381T73 411T99 433T132 442Q161 442 183 430T214 408T225 388Q227 382 228 382T236 389Q284 441 347 441H350Q398 441 422 400Q430 381 430 363Q430 333 417 315T391 292T366 288Q346 288 334 299T322 328Q322 376 378 392Q356 405 342 405Q286 405 239 331Q229 315 224 298T190 165Q156 25 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 114 189T154 366Q154 405 128 405Q107 405 92 377T68 316T57 280Q55 278 41 278H27Q21 284 21 287Z" style="stroke-width:3;"></path></g><g data-mml-node="mn" transform="translate(484,413) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g><g data-mml-node="mo" transform="translate(5308,0)"><path data-c="2B" d="M56 237T56 250T70 270H369V420L370 570Q380 583 389 583Q402 583 409 568V270H707Q722 262 722 250T707 230H409V-68Q401 -82 391 -82H389H387Q375 -82 369 -68V230H70Q56 237 56 250Z" style="stroke-width:3;"></path></g><g data-mml-node="msup" transform="translate(6308.2,0)"><g data-mml-node="mi"><path data-c="1D44E" d="M33 157Q33 258 109 349T280 441Q331 441 370 392Q386 422 416 422Q429 422 439 414T449 394Q449 381 412 234T374 68Q374 43 381 35T402 26Q411 27 422 35Q443 55 463 131Q469 151 473 152Q475 153 483 153H487Q506 153 506 144Q506 138 501 117T481 63T449 13Q436 0 417 -8Q409 -10 393 -10Q359 -10 336 5T306 36L300 51Q299 52 296 50Q294 48 292 46Q233 -10 172 -10Q117 -10 75 30T33 157ZM351 328Q351 334 346 350T323 385T277 405Q242 405 210 374T160 293Q131 214 119 129Q119 126 119 118T118 106Q118 61 136 44T179 26Q217 26 254 59T298 110Q300 114 325 217T351 328Z" style="stroke-width:3;"></path></g><g data-mml-node="mn" transform="translate(562,413) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g><g data-mml-node="msup" transform="translate(7440.4,0)"><g data-mml-node="mi"><path data-c="63" d="M370 305T349 305T313 320T297 358Q297 381 312 396Q317 401 317 402T307 404Q281 408 258 408Q209 408 178 376Q131 329 131 219Q131 137 162 90Q203 29 272 29Q313 29 338 55T374 117Q376 125 379 127T395 129H409Q415 123 415 120Q415 116 411 104T395 71T366 33T318 2T249 -11Q163 -11 99 53T34 214Q34 318 99 383T250 448T370 421T404 357Q404 334 387 320Z" style="stroke-width:3;"></path><path data-c="6F" d="M28 214Q28 309 93 378T250 448Q340 448 405 380T471 215Q471 120 407 55T250 -10Q153 -10 91 57T28 214ZM250 30Q372 30 372 193V225V250Q372 272 371 288T364 326T348 362T317 390T268 410Q263 411 252 411Q222 411 195 399Q152 377 139 338T126 246V226Q126 130 145 91Q177 30 250 30Z" transform="translate(444,0)" style="stroke-width:3;"></path><path data-c="73" d="M295 316Q295 356 268 385T190 414Q154 414 128 401Q98 382 98 349Q97 344 98 336T114 312T157 287Q175 282 201 278T245 269T277 256Q294 248 310 236T342 195T359 133Q359 71 321 31T198 -10H190Q138 -10 94 26L86 19L77 10Q71 4 65 -1L54 -11H46H42Q39 -11 33 -5V74V132Q33 153 35 157T45 162H54Q66 162 70 158T75 146T82 119T101 77Q136 26 198 26Q295 26 295 104Q295 133 277 151Q257 175 194 187T111 210Q75 227 54 256T33 318Q33 357 50 384T93 424T143 442T187 447H198Q238 447 268 432L283 424L292 431Q302 440 314 448H322H326Q329 448 335 442V310L329 304H301Q295 310 295 316Z" transform="translate(944,0)" style="stroke-width:3;"></path></g><g data-mml-node="mn" transform="translate(1371,413) scale(0.707)"><path data-c="32" d="M109 429Q82 429 66 447T50 491Q50 562 103 614T235 666Q326 666 387 610T449 465Q449 422 429 383T381 315T301 241Q265 210 201 149L142 93L218 92Q375 92 385 97Q392 99 409 186V189H449V186Q448 183 436 95T421 3V0H50V19V31Q50 38 56 46T86 81Q115 113 136 137Q145 147 170 174T204 211T233 244T261 278T284 308T305 340T320 369T333 401T340 431T343 464Q343 527 309 573T212 619Q179 619 154 602T119 569T109 550Q109 549 114 549Q132 549 151 535T170 489Q170 464 154 447T109 429Z" style="stroke-width:3;"></path></g></g><g data-mml-node="mo" transform="translate(9215,0)"><path data-c="2061" d="" style="stroke-width:3;"></path></g><g data-mml-node="mi" transform="translate(9381.7,0)"><path data-c="1D703" d="M35 200Q35 302 74 415T180 610T319 704Q320 704 327 704T339 705Q393 701 423 656Q462 596 462 495Q462 380 417 261T302 66T168 -10H161Q125 -10 99 10T60 63T41 130T35 200ZM383 566Q383 668 330 668Q294 668 260 623T204 521T170 421T157 371Q206 370 254 370L351 371Q352 372 359 404T375 484T383 566ZM113 132Q113 26 166 26Q181 26 198 36T239 74T287 161T335 307L340 324H145Q145 321 136 286T120 208T113 132Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(9850.7,0)"><path data-c="2E" d="M78 60Q78 84 95 102T138 120Q162 120 180 104T199 61Q199 36 182 18T139 0T96 17T78 60Z" style="stroke-width:3;"></path></g></g></g>',1)]))),e[21]||(e[21]=t("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[t("mi",{mathvariant:"normal"},"Σ"),t("mo",{stretchy:"false"},"("),t("mi",null,"r"),t("mo",null,","),t("mi",null,"θ"),t("mo",{stretchy:"false"},")"),t("mo",null,"="),t("msup",null,[t("mi",null,"r"),t("mn",null,"2")]),t("mo",null,"+"),t("msup",null,[t("mi",null,"a"),t("mn",null,"2")]),t("msup",null,[t("mi",null,"cos"),t("mn",null,"2")]),t("mo",{"data-mjx-texclass":"NONE"},"⁡"),t("mi",null,"θ"),t("mo",null,".")])],-1))]),e[28]||(e[28]=l('<p><img src="'+o+`" alt=""></p><p style="text-align:center;">Coordinate evolution with Mino time.</p><h3 id="cameras" tabindex="-1">Cameras <a class="header-anchor" href="#cameras" aria-label="Permalink to &quot;Cameras&quot;">​</a></h3><p>Cameras store pre-computed information that is constant for a given camera location. There are currently two types of cameras which can be used for either &#39;slow light&#39; or &#39;fast light&#39; raytracing.</p><ul><li><p><code>IntensityCamera</code> : Precomputes geodesic information necessary to solve the &#39;fast light&#39; raytracing problem.</p></li><li><p><code>SlowLightIntensityCamera</code> : Precomputes geodesic information necessary to solve the &#39;slow light&#39; raytracing problem.</p></li></ul><p>The GPU arrays can be passed to the cameras on construction to raytrace enforce raytracing on the GPU. An sketch of how to do this with a CUDA array is:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CUDA</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> </span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">store </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CUDA</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">fill</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, sze, sze)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">camera </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Krang</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SlowLightIntensityCamera</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(metric, θo, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ρmax, ρmax, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ρmax, ρmax, sze, A</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CuArray)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Krang</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">render!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(store, camera, scene)</span></span></code></pre></div><h2 id="examples" tabindex="-1">Examples <a class="header-anchor" href="#examples" aria-label="Permalink to &quot;Examples&quot;">​</a></h2><p>Some pedagogical example usage of the raytracing code in the <a href="./examples/coordinate-example">Examples</a> section. These examples primarily use <code>CairoMakie</code> and <code>GLMakie</code> for plotting. You can install these by enering the <code>julia</code> package mode and running:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> add CairoMakie GLMakie</span></span></code></pre></div>`,10))])}const F=r(d,[["render",M]]);export{D as __pageData,F as default};