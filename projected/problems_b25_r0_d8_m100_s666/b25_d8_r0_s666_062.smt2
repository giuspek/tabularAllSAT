(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (or A15 A8))) (let ((.def_1 (not .def_0))) (let ((.def_2 (not A7))) (let ((.def_3 (not A0))) (let ((.def_4 (and .def_3 .def_2))) (let ((.def_5 (not .def_4))) (let ((.def_6 (or .def_5 .def_1))) (let ((.def_7 (or A19 A10))) (let ((.def_8 (not .def_7))) (let ((.def_9 (not A19))) (let ((.def_10 (or .def_9 .def_2))) (let ((.def_11 (or .def_10 .def_8))) (let ((.def_12 (or .def_11 .def_6))) (let ((.def_13 (not A4))) (let ((.def_14 (and A8 .def_13))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and A3 A4))) (let ((.def_17 (= .def_16 .def_15))) (let ((.def_18 (not A15))) (let ((.def_19 (not A22))) (let ((.def_20 (and .def_19 .def_18))) (let ((.def_21 (or A21 A10))) (let ((.def_22 (not .def_21))) (let ((.def_23 (and .def_22 .def_20))) (let ((.def_24 (not .def_23))) (let ((.def_25 (or .def_24 .def_17))) (let ((.def_26 (and .def_25 .def_12))) (let ((.def_27 (and A10 A19))) (let ((.def_28 (not A17))) (let ((.def_29 (or A1 .def_28))) (let ((.def_30 (not .def_29))) (let ((.def_31 (or .def_30 .def_27))) (let ((.def_32 (not .def_31))) (let ((.def_33 (not A11))) (let ((.def_34 (or A21 .def_33))) (let ((.def_35 (not .def_34))) (let ((.def_36 (or A2 A9))) (let ((.def_37 (not .def_36))) (let ((.def_38 (or .def_37 .def_35))) (let ((.def_39 (or .def_38 .def_32))) (let ((.def_40 (not .def_39))) (let ((.def_41 (and A23 A2))) (let ((.def_42 (not .def_41))) (let ((.def_43 (not A16))) (let ((.def_44 (or A20 .def_43))) (let ((.def_45 (and .def_44 .def_42))) (let ((.def_46 (not .def_45))) (let ((.def_47 (or .def_28 A23))) (let ((.def_48 (not .def_47))) (let ((.def_49 (not A9))) (let ((.def_50 (not A8))) (let ((.def_51 (and .def_50 .def_49))) (let ((.def_52 (not .def_51))) (let ((.def_53 (or .def_52 .def_48))) (let ((.def_54 (not .def_53))) (let ((.def_55 (and .def_54 .def_46))) (let ((.def_56 (not .def_55))) (let ((.def_57 (and .def_56 .def_40))) (let ((.def_58 (or .def_57 .def_26))) (let ((.def_59 (not A6))) (let ((.def_60 (or .def_59 .def_9))) (let ((.def_61 (not A1))) (let ((.def_62 (= .def_61 .def_3))) (let ((.def_63 (not .def_62))) (let ((.def_64 (= .def_63 .def_60))) (let ((.def_65 (not .def_64))) (let ((.def_66 (and A11 .def_43))) (let ((.def_67 (not A3))) (let ((.def_68 (and .def_3 .def_67))) (let ((.def_69 (and .def_68 .def_66))) (let ((.def_70 (or .def_69 .def_65))) (let ((.def_71 (not .def_70))) (let ((.def_72 (or .def_43 .def_9))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and A18 A13))) (let ((.def_75 (not .def_74))) (let ((.def_76 (or .def_75 .def_73))) (let ((.def_77 (not A12))) (let ((.def_78 (and .def_67 .def_77))) (let ((.def_79 (not .def_78))) (let ((.def_80 (not A20))) (let ((.def_81 (not A14))) (let ((.def_82 (and .def_81 .def_80))) (let ((.def_83 (or .def_82 .def_79))) (let ((.def_84 (not .def_83))) (let ((.def_85 (or .def_84 .def_76))) (let ((.def_86 (not .def_85))) (let ((.def_87 (or .def_86 .def_71))) (let ((.def_88 (and .def_50 .def_19))) (let ((.def_89 (not .def_88))) (let ((.def_90 (and .def_61 .def_13))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and .def_91 .def_89))) (let ((.def_93 (not .def_92))) (let ((.def_94 (and .def_2 .def_50))) (let ((.def_95 (= A22 .def_80))) (let ((.def_96 (not .def_95))) (let ((.def_97 (and .def_96 .def_94))) (let ((.def_98 (not .def_97))) (let ((.def_99 (or .def_98 .def_93))) (let ((.def_100 (= .def_77 .def_3))) (let ((.def_101 (not A23))) (let ((.def_102 (and A15 .def_101))) (let ((.def_103 (not .def_102))) (let ((.def_104 (and .def_103 .def_100))) (let ((.def_105 (and A20 A22))) (let ((.def_106 (not .def_105))) (let ((.def_107 (and .def_101 A5))) (let ((.def_108 (or .def_107 .def_106))) (let ((.def_109 (not .def_108))) (let ((.def_110 (or .def_109 .def_104))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_99))) (let ((.def_113 (not .def_112))) (let ((.def_114 (and .def_113 .def_87))) (let ((.def_115 (not .def_114))) (let ((.def_116 (and .def_115 .def_58))) (let ((.def_117 (and A21 .def_43))) (let ((.def_118 (= .def_18 A2))) (let ((.def_119 (or .def_118 .def_117))) (let ((.def_120 (and .def_49 A24))) (let ((.def_121 (not .def_120))) (let ((.def_122 (or A11 .def_3))) (let ((.def_123 (not .def_122))) (let ((.def_124 (or .def_123 .def_121))) (let ((.def_125 (and .def_124 .def_119))) (let ((.def_126 (not .def_125))) (let ((.def_127 (and A12 A0))) (let ((.def_128 (or A6 A4))) (let ((.def_129 (and .def_128 .def_127))) (let ((.def_130 (not .def_129))) (let ((.def_131 (and A6 A10))) (let ((.def_132 (= .def_80 .def_13))) (let ((.def_133 (not .def_132))) (let ((.def_134 (or .def_133 .def_131))) (let ((.def_135 (not .def_134))) (let ((.def_136 (and .def_135 .def_130))) (let ((.def_137 (not .def_136))) (let ((.def_138 (and .def_137 .def_126))) (let ((.def_139 (not .def_138))) (let ((.def_140 (or A12 A2))) (let ((.def_141 (not .def_140))) (let ((.def_142 (not A5))) (let ((.def_143 (and .def_18 .def_142))) (let ((.def_144 (and .def_143 .def_141))) (let ((.def_145 (or A8 A16))) (let ((.def_146 (not .def_145))) (let ((.def_147 (not A2))) (let ((.def_148 (and .def_147 .def_61))) (let ((.def_149 (not .def_148))) (let ((.def_150 (or .def_149 .def_146))) (let ((.def_151 (= .def_150 .def_144))) (let ((.def_152 (= .def_43 A13))) (let ((.def_153 (not .def_152))) (let ((.def_154 (and .def_2 .def_43))) (let ((.def_155 (not .def_154))) (let ((.def_156 (= .def_155 .def_153))) (let ((.def_157 (or A24 A18))) (let ((.def_158 (not .def_157))) (let ((.def_159 (and .def_101 A3))) (let ((.def_160 (or .def_159 .def_158))) (let ((.def_161 (not .def_160))) (let ((.def_162 (= .def_161 .def_156))) (let ((.def_163 (not .def_162))) (let ((.def_164 (or .def_163 .def_151))) (let ((.def_165 (or .def_164 .def_139))) (let ((.def_166 (not A13))) (let ((.def_167 (and A19 .def_166))) (let ((.def_168 (not .def_167))) (let ((.def_169 (or A8 .def_61))) (let ((.def_170 (and .def_169 .def_168))) (let ((.def_171 (and A23 .def_3))) (let ((.def_172 (and .def_50 A13))) (let ((.def_173 (not .def_172))) (let ((.def_174 (and .def_173 .def_171))) (let ((.def_175 (or .def_174 .def_170))) (let ((.def_176 (not .def_175))) (let ((.def_177 (and A8 .def_33))) (let ((.def_178 (and .def_77 A9))) (let ((.def_179 (or .def_178 .def_177))) (let ((.def_180 (and A9 .def_80))) (let ((.def_181 (not .def_180))) (let ((.def_182 (and A22 A7))) (let ((.def_183 (not .def_182))) (let ((.def_184 (and .def_183 .def_181))) (let ((.def_185 (or .def_184 .def_179))) (let ((.def_186 (not .def_185))) (let ((.def_187 (and .def_186 .def_176))) (let ((.def_188 (not .def_187))) (let ((.def_189 (and A8 A4))) (let ((.def_190 (not .def_189))) (let ((.def_191 (not A18))) (let ((.def_192 (not A24))) (let ((.def_193 (or .def_192 .def_191))) (let ((.def_194 (or .def_193 .def_190))) (let ((.def_195 (and .def_18 .def_192))) (let ((.def_196 (not .def_195))) (let ((.def_197 (= A0 .def_77))) (let ((.def_198 (not .def_197))) (let ((.def_199 (and .def_198 .def_196))) (let ((.def_200 (and .def_199 .def_194))) (let ((.def_201 (not .def_200))) (let ((.def_202 (and .def_3 .def_77))) (let ((.def_203 (or .def_166 A16))) (let ((.def_204 (and .def_203 .def_202))) (let ((.def_205 (or .def_61 .def_166))) (let ((.def_206 (not .def_205))) (let ((.def_207 (= .def_9 A17))) (let ((.def_208 (and .def_207 .def_206))) (let ((.def_209 (not .def_208))) (let ((.def_210 (or .def_209 .def_204))) (let ((.def_211 (not .def_210))) (let ((.def_212 (and .def_211 .def_201))) (let ((.def_213 (and .def_212 .def_188))) (let ((.def_214 (and .def_213 .def_165))) (let ((.def_215 (or .def_214 .def_116))) (let ((.def_216 (and .def_142 .def_2))) (let ((.def_217 (or .def_50 .def_101))) (let ((.def_218 (not .def_217))) (let ((.def_219 (and .def_218 .def_216))) (let ((.def_220 (not .def_219))) (let ((.def_221 (or A24 .def_77))) (let ((.def_222 (not .def_221))) (let ((.def_223 (and A2 A1))) (let ((.def_224 (not .def_223))) (let ((.def_225 (and .def_224 .def_222))) (let ((.def_226 (not .def_225))) (let ((.def_227 (and .def_226 .def_220))) (let ((.def_228 (or A13 A9))) (let ((.def_229 (= .def_49 .def_77))) (let ((.def_230 (and .def_229 .def_228))) (let ((.def_231 (not .def_230))) (let ((.def_232 (not A10))) (let ((.def_233 (and A14 .def_232))) (let ((.def_234 (not .def_233))) (let ((.def_235 (or A1 A3))) (let ((.def_236 (and .def_235 .def_234))) (let ((.def_237 (not .def_236))) (let ((.def_238 (and .def_237 .def_231))) (let ((.def_239 (not .def_238))) (let ((.def_240 (or .def_239 .def_227))) (let ((.def_241 (and .def_67 .def_81))) (let ((.def_242 (not .def_216))) (let ((.def_243 (or .def_242 .def_241))) (let ((.def_244 (or A3 .def_19))) (let ((.def_245 (not .def_244))) (let ((.def_246 (or A8 A1))) (let ((.def_247 (or .def_246 .def_245))) (let ((.def_248 (or .def_247 .def_243))) (let ((.def_249 (and .def_101 A13))) (let ((.def_250 (not .def_249))) (let ((.def_251 (or .def_101 .def_19))) (let ((.def_252 (not .def_251))) (let ((.def_253 (and .def_252 .def_250))) (let ((.def_254 (not .def_253))) (let ((.def_255 (= .def_67 A1))) (let ((.def_256 (and .def_255 .def_78))) (let ((.def_257 (not .def_256))) (let ((.def_258 (and .def_257 .def_254))) (let ((.def_259 (and .def_258 .def_248))) (let ((.def_260 (not .def_259))) (let ((.def_261 (and .def_260 .def_240))) (let ((.def_262 (not .def_261))) (let ((.def_263 (and A21 A16))) (let ((.def_264 (or A16 .def_77))) (let ((.def_265 (not .def_264))) (let ((.def_266 (or .def_265 .def_263))) (let ((.def_267 (not .def_266))) (let ((.def_268 (= .def_43 .def_3))) (let ((.def_269 (and .def_61 .def_19))) (let ((.def_270 (not .def_269))) (let ((.def_271 (or .def_270 .def_268))) (let ((.def_272 (not .def_271))) (let ((.def_273 (or .def_272 .def_267))) (let ((.def_274 (or .def_191 .def_77))) (let ((.def_275 (not .def_274))) (let ((.def_276 (and .def_80 .def_18))) (let ((.def_277 (not .def_276))) (let ((.def_278 (or .def_277 .def_275))) (let ((.def_279 (or A9 A21))) (let ((.def_280 (not .def_279))) (let ((.def_281 (= A0 A14))) (let ((.def_282 (not .def_281))) (let ((.def_283 (or .def_282 .def_280))) (let ((.def_284 (not .def_283))) (let ((.def_285 (or .def_284 .def_278))) (let ((.def_286 (and .def_285 .def_273))) (let ((.def_287 (and .def_43 A24))) (let ((.def_288 (= .def_107 .def_287))) (let ((.def_289 (not .def_288))) (let ((.def_290 (and .def_192 A22))) (let ((.def_291 (not .def_290))) (let ((.def_292 (or .def_77 A9))) (let ((.def_293 (and .def_292 .def_291))) (let ((.def_294 (or .def_293 .def_289))) (let ((.def_295 (or A8 A17))) (let ((.def_296 (not .def_295))) (let ((.def_297 (= A5 A7))) (let ((.def_298 (not .def_297))) (let ((.def_299 (and .def_298 .def_296))) (let ((.def_300 (not .def_299))) (let ((.def_301 (or .def_9 A18))) (let ((.def_302 (or .def_3 A18))) (let ((.def_303 (and .def_302 .def_301))) (let ((.def_304 (and .def_303 .def_300))) (let ((.def_305 (not .def_304))) (let ((.def_306 (or .def_305 .def_294))) (let ((.def_307 (= .def_306 .def_286))) (let ((.def_308 (or .def_307 .def_262))) (let ((.def_309 (and A9 .def_67))) (let ((.def_310 (or .def_33 A6))) (let ((.def_311 (not .def_310))) (let ((.def_312 (or .def_311 .def_309))) (let ((.def_313 (not .def_312))) (let ((.def_314 (and A11 .def_13))) (let ((.def_315 (not .def_314))) (let ((.def_316 (not A21))) (let ((.def_317 (and .def_316 A3))) (let ((.def_318 (and .def_317 .def_315))) (let ((.def_319 (not .def_318))) (let ((.def_320 (= .def_319 .def_313))) (let ((.def_321 (or .def_9 A10))) (let ((.def_322 (and A17 .def_2))) (let ((.def_323 (not .def_322))) (let ((.def_324 (or .def_323 .def_321))) (let ((.def_325 (and .def_77 .def_59))) (let ((.def_326 (not .def_325))) (let ((.def_327 (or A15 .def_33))) (let ((.def_328 (not .def_327))) (let ((.def_329 (= .def_328 .def_326))) (let ((.def_330 (not .def_329))) (let ((.def_331 (and .def_330 .def_324))) (let ((.def_332 (and .def_331 .def_320))) (let ((.def_333 (not .def_332))) (let ((.def_334 (and .def_80 A11))) (let ((.def_335 (or A6 A0))) (let ((.def_336 (not .def_335))) (let ((.def_337 (or .def_336 .def_334))) (let ((.def_338 (and A14 .def_19))) (let ((.def_339 (not .def_338))) (let ((.def_340 (and .def_101 A23))) (let ((.def_341 (not .def_340))) (let ((.def_342 (and .def_341 .def_339))) (let ((.def_343 (not .def_342))) (let ((.def_344 (or .def_343 .def_337))) (let ((.def_345 (not .def_344))) (let ((.def_346 (and .def_147 .def_192))) (let ((.def_347 (not .def_346))) (let ((.def_348 (or A4 A1))) (let ((.def_349 (not .def_348))) (let ((.def_350 (= .def_349 .def_347))) (let ((.def_351 (not .def_350))) (let ((.def_352 (and A1 .def_50))) (let ((.def_353 (or .def_2 A11))) (let ((.def_354 (not .def_353))) (let ((.def_355 (or .def_354 .def_352))) (let ((.def_356 (not .def_355))) (let ((.def_357 (or .def_356 .def_351))) (let ((.def_358 (and .def_357 .def_345))) (let ((.def_359 (= .def_358 .def_333))) (let ((.def_360 (not .def_359))) (let ((.def_361 (or .def_80 A5))) (let ((.def_362 (not .def_361))) (let ((.def_363 (and A14 A7))) (let ((.def_364 (not .def_363))) (let ((.def_365 (and .def_364 .def_362))) (let ((.def_366 (= A7 .def_80))) (let ((.def_367 (or A21 A4))) (let ((.def_368 (not .def_367))) (let ((.def_369 (or .def_368 .def_366))) (let ((.def_370 (not .def_369))) (let ((.def_371 (or .def_370 .def_365))) (let ((.def_372 (= A8 A14))) (let ((.def_373 (not .def_372))) (let ((.def_374 (or .def_59 .def_2))) (let ((.def_375 (and .def_374 .def_373))) (let ((.def_376 (not .def_375))) (let ((.def_377 (and .def_191 .def_191))) (let ((.def_378 (not .def_377))) (let ((.def_379 (or A5 A21))) (let ((.def_380 (or .def_379 .def_378))) (let ((.def_381 (not .def_380))) (let ((.def_382 (and .def_381 .def_376))) (let ((.def_383 (not .def_382))) (let ((.def_384 (and .def_383 .def_371))) (let ((.def_385 (or .def_192 .def_18))) (let ((.def_386 (not .def_385))) (let ((.def_387 (= A18 A22))) (let ((.def_388 (not .def_387))) (let ((.def_389 (and .def_388 .def_386))) (let ((.def_390 (or A3 .def_81))) (let ((.def_391 (or A12 .def_192))) (let ((.def_392 (not .def_391))) (let ((.def_393 (and .def_392 .def_390))) (let ((.def_394 (or .def_393 .def_389))) (let ((.def_395 (or A24 A12))) (let ((.def_396 (not .def_395))) (let ((.def_397 (and A2 .def_80))) (let ((.def_398 (or .def_397 .def_396))) (let ((.def_399 (and .def_28 .def_28))) (let ((.def_400 (not .def_399))) (let ((.def_401 (and A21 .def_13))) (let ((.def_402 (and .def_401 .def_400))) (let ((.def_403 (not .def_402))) (let ((.def_404 (or .def_403 .def_398))) (let ((.def_405 (not .def_404))) (let ((.def_406 (or .def_405 .def_394))) (let ((.def_407 (not .def_406))) (let ((.def_408 (and .def_407 .def_384))) (let ((.def_409 (not .def_408))) (let ((.def_410 (or .def_409 .def_360))) (let ((.def_411 (or .def_410 .def_308))) (let ((.def_412 (or .def_411 .def_215))) .def_412))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)