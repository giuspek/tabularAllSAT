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
(assert (let ((.def_0 (not A5))) (let ((.def_1 (or A5 .def_0))) (let ((.def_2 (not A10))) (let ((.def_3 (or .def_2 A11))) (let ((.def_4 (and .def_3 .def_1))) (let ((.def_5 (not .def_4))) (let ((.def_6 (not A19))) (let ((.def_7 (not A1))) (let ((.def_8 (or .def_7 .def_6))) (let ((.def_9 (and A19 A11))) (let ((.def_10 (not .def_9))) (let ((.def_11 (or .def_10 .def_8))) (let ((.def_12 (and .def_11 .def_5))) (let ((.def_13 (not .def_12))) (let ((.def_14 (not A4))) (let ((.def_15 (not A11))) (let ((.def_16 (and .def_15 .def_14))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A23))) (let ((.def_19 (or A1 .def_18))) (let ((.def_20 (not .def_19))) (let ((.def_21 (and .def_20 .def_17))) (let ((.def_22 (not .def_21))) (let ((.def_23 (or A16 A23))) (let ((.def_24 (not .def_23))) (let ((.def_25 (and A8 A5))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and .def_26 .def_24))) (let ((.def_28 (and .def_27 .def_22))) (let ((.def_29 (= .def_28 .def_13))) (let ((.def_30 (or A5 A14))) (let ((.def_31 (or A14 A5))) (let ((.def_32 (or .def_31 .def_30))) (let ((.def_33 (not .def_32))) (let ((.def_34 (not A18))) (let ((.def_35 (or A14 .def_34))) (let ((.def_36 (not .def_35))) (let ((.def_37 (not A16))) (let ((.def_38 (or A13 .def_37))) (let ((.def_39 (not .def_38))) (let ((.def_40 (and .def_39 .def_36))) (let ((.def_41 (or .def_40 .def_33))) (let ((.def_42 (not .def_41))) (let ((.def_43 (not A2))) (let ((.def_44 (not A9))) (let ((.def_45 (or .def_44 .def_43))) (let ((.def_46 (not .def_45))) (let ((.def_47 (not A21))) (let ((.def_48 (or A6 .def_47))) (let ((.def_49 (not .def_48))) (let ((.def_50 (and .def_49 .def_46))) (let ((.def_51 (not .def_50))) (let ((.def_52 (not A13))) (let ((.def_53 (or .def_52 .def_2))) (let ((.def_54 (not A3))) (let ((.def_55 (or .def_54 A24))) (let ((.def_56 (not .def_55))) (let ((.def_57 (and .def_56 .def_53))) (let ((.def_58 (and .def_57 .def_51))) (let ((.def_59 (not .def_58))) (let ((.def_60 (and .def_59 .def_42))) (let ((.def_61 (not .def_60))) (let ((.def_62 (or .def_61 .def_29))) (let ((.def_63 (not .def_62))) (let ((.def_64 (or .def_54 A1))) (let ((.def_65 (not .def_64))) (let ((.def_66 (and .def_54 A7))) (let ((.def_67 (or .def_66 .def_65))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_52 A9))) (let ((.def_70 (not A15))) (let ((.def_71 (or A8 .def_70))) (let ((.def_72 (not .def_71))) (let ((.def_73 (= .def_72 .def_69))) (let ((.def_74 (= .def_73 .def_68))) (let ((.def_75 (and .def_43 A5))) (let ((.def_76 (not A7))) (let ((.def_77 (= .def_76 .def_2))) (let ((.def_78 (or .def_77 .def_75))) (let ((.def_79 (and .def_2 A14))) (let ((.def_80 (or .def_34 A17))) (let ((.def_81 (not .def_80))) (let ((.def_82 (or .def_81 .def_79))) (let ((.def_83 (not .def_82))) (let ((.def_84 (and .def_83 .def_78))) (let ((.def_85 (= .def_84 .def_74))) (let ((.def_86 (or .def_52 .def_70))) (let ((.def_87 (not .def_86))) (let ((.def_88 (or A20 .def_7))) (let ((.def_89 (not .def_88))) (let ((.def_90 (= .def_89 .def_87))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and .def_54 .def_14))) (let ((.def_93 (not A6))) (let ((.def_94 (or A6 .def_93))) (let ((.def_95 (not .def_94))) (let ((.def_96 (or .def_95 .def_92))) (let ((.def_97 (or .def_96 .def_91))) (let ((.def_98 (not A20))) (let ((.def_99 (and .def_76 .def_98))) (let ((.def_100 (and .def_70 A19))) (let ((.def_101 (or .def_100 .def_99))) (let ((.def_102 (not .def_101))) (let ((.def_103 (or A6 A8))) (let ((.def_104 (not .def_103))) (let ((.def_105 (and A21 A16))) (let ((.def_106 (or .def_105 .def_104))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_107 .def_102))) (let ((.def_109 (not .def_108))) (let ((.def_110 (or .def_109 .def_97))) (let ((.def_111 (not .def_110))) (let ((.def_112 (and .def_111 .def_85))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_113 .def_63))) (let ((.def_115 (and .def_76 A22))) (let ((.def_116 (not .def_115))) (let ((.def_117 (or .def_44 .def_70))) (let ((.def_118 (not .def_117))) (let ((.def_119 (and .def_118 .def_116))) (let ((.def_120 (= .def_14 A18))) (let ((.def_121 (or .def_52 .def_15))) (let ((.def_122 (and .def_121 .def_120))) (let ((.def_123 (not .def_122))) (let ((.def_124 (and .def_123 .def_119))) (let ((.def_125 (or A3 .def_47))) (let ((.def_126 (= .def_76 A19))) (let ((.def_127 (or .def_126 .def_125))) (let ((.def_128 (and .def_34 A8))) (let ((.def_129 (not .def_128))) (let ((.def_130 (and A23 A15))) (let ((.def_131 (and .def_130 .def_129))) (let ((.def_132 (not .def_131))) (let ((.def_133 (or .def_132 .def_127))) (let ((.def_134 (not .def_133))) (let ((.def_135 (and .def_134 .def_124))) (let ((.def_136 (not A17))) (let ((.def_137 (and .def_136 .def_44))) (let ((.def_138 (not A22))) (let ((.def_139 (or .def_138 .def_98))) (let ((.def_140 (not .def_139))) (let ((.def_141 (and .def_140 .def_137))) (let ((.def_142 (not .def_141))) (let ((.def_143 (and .def_34 .def_18))) (let ((.def_144 (not .def_143))) (let ((.def_145 (not A0))) (let ((.def_146 (and A15 .def_145))) (let ((.def_147 (not .def_146))) (let ((.def_148 (and .def_147 .def_144))) (let ((.def_149 (not .def_148))) (let ((.def_150 (and .def_149 .def_142))) (let ((.def_151 (not .def_150))) (let ((.def_152 (and .def_145 A20))) (let ((.def_153 (not A14))) (let ((.def_154 (or A2 .def_153))) (let ((.def_155 (or .def_154 .def_152))) (let ((.def_156 (not .def_155))) (let ((.def_157 (and A9 A21))) (let ((.def_158 (not .def_157))) (let ((.def_159 (and .def_43 A20))) (let ((.def_160 (and .def_159 .def_158))) (let ((.def_161 (= .def_160 .def_156))) (let ((.def_162 (not .def_161))) (let ((.def_163 (and .def_162 .def_151))) (let ((.def_164 (or .def_163 .def_135))) (let ((.def_165 (not .def_164))) (let ((.def_166 (or .def_44 A5))) (let ((.def_167 (not .def_166))) (let ((.def_168 (= .def_37 A14))) (let ((.def_169 (or .def_168 .def_167))) (let ((.def_170 (not .def_169))) (let ((.def_171 (= A0 A2))) (let ((.def_172 (not .def_171))) (let ((.def_173 (and .def_7 A9))) (let ((.def_174 (not .def_173))) (let ((.def_175 (or .def_174 .def_172))) (let ((.def_176 (and .def_175 .def_170))) (let ((.def_177 (not .def_176))) (let ((.def_178 (and .def_145 .def_18))) (let ((.def_179 (not .def_178))) (let ((.def_180 (or .def_6 A16))) (let ((.def_181 (not .def_180))) (let ((.def_182 (or .def_181 .def_179))) (let ((.def_183 (or .def_7 .def_47))) (let ((.def_184 (not .def_183))) (let ((.def_185 (or .def_153 .def_44))) (let ((.def_186 (not .def_185))) (let ((.def_187 (or .def_186 .def_184))) (let ((.def_188 (not .def_187))) (let ((.def_189 (or .def_188 .def_182))) (let ((.def_190 (not .def_189))) (let ((.def_191 (or .def_190 .def_177))) (let ((.def_192 (not .def_191))) (let ((.def_193 (and .def_93 .def_138))) (let ((.def_194 (or .def_47 A24))) (let ((.def_195 (not .def_194))) (let ((.def_196 (or .def_195 .def_193))) (let ((.def_197 (not .def_196))) (let ((.def_198 (or .def_34 A4))) (let ((.def_199 (not .def_198))) (let ((.def_200 (= .def_145 .def_98))) (let ((.def_201 (or .def_200 .def_199))) (let ((.def_202 (or .def_201 .def_197))) (let ((.def_203 (and .def_153 A12))) (let ((.def_204 (and A3 .def_34))) (let ((.def_205 (not .def_204))) (let ((.def_206 (or .def_205 .def_203))) (let ((.def_207 (not A8))) (let ((.def_208 (or A1 .def_207))) (let ((.def_209 (and A16 A5))) (let ((.def_210 (not .def_209))) (let ((.def_211 (and .def_210 .def_208))) (let ((.def_212 (and .def_211 .def_206))) (let ((.def_213 (not .def_212))) (let ((.def_214 (and .def_213 .def_202))) (let ((.def_215 (not .def_214))) (let ((.def_216 (or .def_215 .def_192))) (let ((.def_217 (and .def_216 .def_165))) (let ((.def_218 (not .def_217))) (let ((.def_219 (and .def_218 .def_114))) (let ((.def_220 (not A12))) (let ((.def_221 (or .def_76 .def_220))) (let ((.def_222 (and A13 A23))) (let ((.def_223 (not .def_222))) (let ((.def_224 (and .def_223 .def_221))) (let ((.def_225 (or .def_76 .def_2))) (let ((.def_226 (not .def_225))) (let ((.def_227 (and A0 A3))) (let ((.def_228 (and .def_227 .def_226))) (let ((.def_229 (and .def_228 .def_224))) (let ((.def_230 (and .def_15 .def_54))) (let ((.def_231 (not .def_230))) (let ((.def_232 (and .def_52 A6))) (let ((.def_233 (or .def_232 .def_231))) (let ((.def_234 (or .def_138 A17))) (let ((.def_235 (or A8 A18))) (let ((.def_236 (not .def_235))) (let ((.def_237 (or .def_236 .def_234))) (let ((.def_238 (not .def_237))) (let ((.def_239 (and .def_238 .def_233))) (let ((.def_240 (not .def_239))) (let ((.def_241 (= .def_240 .def_229))) (let ((.def_242 (and A19 .def_153))) (let ((.def_243 (not .def_242))) (let ((.def_244 (and .def_43 .def_220))) (let ((.def_245 (not .def_244))) (let ((.def_246 (and .def_245 .def_243))) (let ((.def_247 (not .def_246))) (let ((.def_248 (and A5 .def_145))) (let ((.def_249 (and .def_18 .def_98))) (let ((.def_250 (not .def_249))) (let ((.def_251 (and .def_250 .def_248))) (let ((.def_252 (not .def_251))) (let ((.def_253 (or .def_252 .def_247))) (let ((.def_254 (not .def_253))) (let ((.def_255 (= .def_145 A5))) (let ((.def_256 (or .def_98 A9))) (let ((.def_257 (not .def_256))) (let ((.def_258 (and .def_257 .def_255))) (let ((.def_259 (or .def_145 A24))) (let ((.def_260 (not .def_259))) (let ((.def_261 (and .def_43 .def_138))) (let ((.def_262 (not .def_261))) (let ((.def_263 (or .def_262 .def_260))) (let ((.def_264 (or .def_263 .def_258))) (let ((.def_265 (not .def_264))) (let ((.def_266 (= .def_265 .def_254))) (let ((.def_267 (not .def_266))) (let ((.def_268 (or .def_267 .def_241))) (let ((.def_269 (and .def_34 A1))) (let ((.def_270 (= .def_138 .def_54))) (let ((.def_271 (not .def_270))) (let ((.def_272 (= .def_271 .def_269))) (let ((.def_273 (or .def_15 .def_7))) (let ((.def_274 (or .def_15 .def_220))) (let ((.def_275 (not .def_274))) (let ((.def_276 (or .def_275 .def_273))) (let ((.def_277 (and .def_276 .def_272))) (let ((.def_278 (not .def_277))) (let ((.def_279 (or A17 A22))) (let ((.def_280 (not .def_279))) (let ((.def_281 (or A19 .def_14))) (let ((.def_282 (not .def_281))) (let ((.def_283 (and .def_282 .def_280))) (let ((.def_284 (not .def_283))) (let ((.def_285 (and .def_138 .def_138))) (let ((.def_286 (or .def_43 .def_145))) (let ((.def_287 (= .def_286 .def_285))) (let ((.def_288 (= .def_287 .def_284))) (let ((.def_289 (not .def_288))) (let ((.def_290 (and .def_289 .def_278))) (let ((.def_291 (not A24))) (let ((.def_292 (and A22 .def_291))) (let ((.def_293 (and A11 A19))) (let ((.def_294 (or .def_293 .def_292))) (let ((.def_295 (and A22 .def_47))) (let ((.def_296 (not .def_295))) (let ((.def_297 (and .def_52 A13))) (let ((.def_298 (not .def_297))) (let ((.def_299 (or .def_298 .def_296))) (let ((.def_300 (and .def_299 .def_294))) (let ((.def_301 (or .def_44 .def_2))) (let ((.def_302 (and A20 .def_145))) (let ((.def_303 (and .def_302 .def_301))) (let ((.def_304 (not .def_303))) (let ((.def_305 (and .def_47 A16))) (let ((.def_306 (not .def_305))) (let ((.def_307 (= .def_207 A0))) (let ((.def_308 (or .def_307 .def_306))) (let ((.def_309 (and .def_308 .def_304))) (let ((.def_310 (not .def_309))) (let ((.def_311 (and .def_310 .def_300))) (let ((.def_312 (not .def_311))) (let ((.def_313 (or .def_312 .def_290))) (let ((.def_314 (not .def_313))) (let ((.def_315 (or .def_314 .def_268))) (let ((.def_316 (not .def_315))) (let ((.def_317 (and .def_7 .def_7))) (let ((.def_318 (not .def_317))) (let ((.def_319 (and A4 A22))) (let ((.def_320 (and .def_319 .def_318))) (let ((.def_321 (not .def_320))) (let ((.def_322 (and A8 .def_52))) (let ((.def_323 (= .def_52 A18))) (let ((.def_324 (and .def_323 .def_322))) (let ((.def_325 (not .def_324))) (let ((.def_326 (and .def_325 .def_321))) (let ((.def_327 (and A24 A9))) (let ((.def_328 (not .def_327))) (let ((.def_329 (= .def_6 A6))) (let ((.def_330 (not .def_329))) (let ((.def_331 (and .def_330 .def_328))) (let ((.def_332 (not .def_331))) (let ((.def_333 (and A5 .def_2))) (let ((.def_334 (not .def_333))) (let ((.def_335 (or .def_54 .def_138))) (let ((.def_336 (not .def_335))) (let ((.def_337 (or .def_336 .def_334))) (let ((.def_338 (and .def_337 .def_332))) (let ((.def_339 (not .def_338))) (let ((.def_340 (= .def_339 .def_326))) (let ((.def_341 (not .def_340))) (let ((.def_342 (or .def_2 .def_70))) (let ((.def_343 (not .def_342))) (let ((.def_344 (= A22 .def_136))) (let ((.def_345 (not .def_344))) (let ((.def_346 (and .def_345 .def_343))) (let ((.def_347 (and .def_7 A24))) (let ((.def_348 (not .def_347))) (let ((.def_349 (= A22 .def_43))) (let ((.def_350 (or .def_349 .def_348))) (let ((.def_351 (not .def_350))) (let ((.def_352 (or .def_351 .def_346))) (let ((.def_353 (= A16 A14))) (let ((.def_354 (not .def_353))) (let ((.def_355 (and A5 .def_47))) (let ((.def_356 (not .def_355))) (let ((.def_357 (or .def_356 .def_354))) (let ((.def_358 (or A10 .def_138))) (let ((.def_359 (not .def_358))) (let ((.def_360 (and A11 .def_207))) (let ((.def_361 (not .def_360))) (let ((.def_362 (and .def_361 .def_359))) (let ((.def_363 (not .def_362))) (let ((.def_364 (= .def_363 .def_357))) (let ((.def_365 (and .def_364 .def_352))) (let ((.def_366 (or .def_365 .def_341))) (let ((.def_367 (and .def_220 A20))) (let ((.def_368 (not .def_367))) (let ((.def_369 (and A5 .def_207))) (let ((.def_370 (or .def_369 .def_368))) (let ((.def_371 (not .def_370))) (let ((.def_372 (and A7 .def_7))) (let ((.def_373 (or A5 A11))) (let ((.def_374 (not .def_373))) (let ((.def_375 (or .def_374 .def_372))) (let ((.def_376 (or .def_375 .def_371))) (let ((.def_377 (not .def_376))) (let ((.def_378 (or A15 .def_6))) (let ((.def_379 (and A14 .def_136))) (let ((.def_380 (not .def_379))) (let ((.def_381 (and .def_380 .def_378))) (let ((.def_382 (or .def_291 .def_44))) (let ((.def_383 (not .def_382))) (let ((.def_384 (= .def_18 A19))) (let ((.def_385 (= .def_384 .def_383))) (let ((.def_386 (not .def_385))) (let ((.def_387 (or .def_386 .def_381))) (let ((.def_388 (and .def_387 .def_377))) (let ((.def_389 (not .def_388))) (let ((.def_390 (= A5 .def_220))) (let ((.def_391 (and .def_291 A17))) (let ((.def_392 (= .def_391 .def_390))) (let ((.def_393 (and .def_47 .def_98))) (let ((.def_394 (and .def_220 A1))) (let ((.def_395 (or .def_394 .def_393))) (let ((.def_396 (not .def_395))) (let ((.def_397 (or .def_396 .def_392))) (let ((.def_398 (not .def_397))) (let ((.def_399 (or A10 .def_37))) (let ((.def_400 (or .def_153 A24))) (let ((.def_401 (not .def_400))) (let ((.def_402 (or .def_401 .def_399))) (let ((.def_403 (not .def_402))) (let ((.def_404 (and .def_43 .def_0))) (let ((.def_405 (not .def_404))) (let ((.def_406 (= .def_138 .def_37))) (let ((.def_407 (not .def_406))) (let ((.def_408 (and .def_407 .def_405))) (let ((.def_409 (not .def_408))) (let ((.def_410 (or .def_409 .def_403))) (let ((.def_411 (or .def_410 .def_398))) (let ((.def_412 (not .def_411))) (let ((.def_413 (or .def_412 .def_389))) (let ((.def_414 (and .def_413 .def_366))) (let ((.def_415 (or .def_414 .def_316))) (let ((.def_416 (not .def_415))) (let ((.def_417 (and .def_416 .def_219))) (let ((.def_418 (not .def_417))) .def_418))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)