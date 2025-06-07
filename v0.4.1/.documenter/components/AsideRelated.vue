
<template>
  <br>
  <h2> Related Software </h2>
  <a
    class="enjoyer"
    href="https://github.com/EHTJulia"
    target="_blank"
  >
    <img width="32" height="32" src="https://avatars.githubusercontent.com/u/102935150?s=200&v=4" />
    <span>
      <p class="extra-info">EHT Julia</p>
      <p class="heading">EHT Julia Organization</p>
      <p class="extra-info">Event Horizon Telescope (EHT) Julia software ecosystem</p>
    </span>
  </a>
  <a
    class="enjoyer"
    href="https://ptiede.github.io/Comrade.jl/stable/"
    target="_blank"
  >
    <span>
      <p class="extra-info">Bayesian inference for VLBI</p>
      <p class="heading">Comrade.jl</p>
      <p class="extra-info">Composable Modeling of Radio Emission</p>
    </span>
  </a>
</template>

<style>
.enjoyer {
  margin-top: 0.5rem;
  margin-bottom: 0.0rem;
  border-radius: 14px;
  padding-top: 0.2rem;
  padding-bottom: 0.2rem;
  position: relative;
  font-size: 0.9rem;
  font-weight: 700;
  line-height: 1.1rem;
  display: flex;
  align-items: center;
  justify-content: center;
  width: 100%;
  gap: 1rem;
  background-color: var(--vp-c-bg-alt);
  border: 2px solid var(--vp-c-bg-alt);
  transition: border-color 0.5s;
}
.enjoyer:hover {
  border: 2px solid var(--vp-c-brand-light);
}
.enjoyer img {
  transition: transform 0.5s;
  transform: scale(1.25);
}
.enjoyer:hover img {
  transform: scale(1.75);
}
.enjoyer .heading {
  background-image: linear-gradient(
    color-mix(in srgb, var(--krang-orange)90%, white),
    color-mix(in srgb, var(--krang-orange)90%, white)
  );
  background-clip: text;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
}
.dark .enjoyer .heading {
  background-image: linear-gradient(
    var(--vp-c-brand-light),
    var(--vp-c-brand-light)
  );
  background-clip: text;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
}
.enjoyer .extra-info {
  color: var(--vp-c-text-1);
  opacity: 0;
  font-size: 0.7rem;
  padding-left: 0.1rem;
  transition: opacity 0.5s;
}
.enjoyer:hover .extra-info {
  opacity: 0.9;
}
</style>
